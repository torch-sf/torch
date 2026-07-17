# ------- GENERATE POLARIS SIMULATION OF A TORCH SNAPSHOT ------- #
# First Written by Brooke Polak, January 2025
# Edited by Donglin Wu, July 2026

import os
import numpy as np
import yt
from amuse.lab import SeBa, Particles

from amuse.units import units
from spherical_detectors import save_angle_detectors

from amuse import io

from scipy.interpolate import RegularGridInterpolator


import warnings


L_sun = 3.828e+26
pc_to_AU = 206264.806




class Polaris:

	# Class constructor

	# sim_directory: path to the POLARIS output directory
	# gridfile: path to the grid file
	# n_bins: number of bins used for midplane files, irrelevant for most applications
	# conv_dens: density conversion factor, set to 1e6 if the grid is in cgs
	# conv_len: length conversion factor, set to 0.01 if the grid is in cgs
	# mass_fraction: set to None to use the dust densities in the grid
	# 				 set to a float if the grid does not have dust densities, and POLARIS will use this value as the dust-to-gas ratio
	# 				 mostly deprecated, as the grid should have dust densities
	# num_threads: number of processors to use for POLARIS, set to -1 to use all available processors
	# subpixel_lvl: maximum subpixel level for POLARIS, set to 3 by default
	# 				for infrared, it is recommended to increase this and do a convergence test, as the subpixel level can affect the accuracy of the results
	
	def __init__(self, sim_directory, gridfile, n_bins="256", conv_dens="1e6",
				 conv_len="0.01", mu="1.3", mass_fraction=None, num_threads="-1",
				 subpixel_lvl="3", 
				 ):

		self.directory = sim_directory
		# make output directories
		if not os.path.isdir(self.directory):
			os.mkdir(self.directory)

		self.gridfile = gridfile
		self.n_bins = n_bins
		self.conv_dens = conv_dens
		self.conv_len = conv_len
		self.mu = mu
		self.mass_fraction = mass_fraction
		self.num_threads = num_threads
		self.dust_components = []
		self.detectors = []
		self.nside_spherical = None
		# create chunk of sources string once to avoid multiple for loops
		self.sources_string = ""

		self.subpixel_lvl = subpixel_lvl

		warnings.warn("The level of subpixel_lvl is set to "
					  +str(self.subpixel_lvl)
					  +". For infrared, it is recommended to do a convergence test on this parameter, as the subpixel level can affect the accuracy of the results.")

		return 


	# Npixel is the number of pixel, λmin is the shortest wavelength, 
	# λmax is the longest wavelength,Nλ is the amount of wavelengths 
	# used for simulation (if 1, only λmin is used), α_x and α_y are 
	# the rotation angles around the first and second rotation axis, 
	# and D is the distance to the observer in meters.
    # dx, dy: side length of the region that the detector looks at (default: the entire grid)
    # delta_x, delta_y: shift of the detector from the center (default: the detector points at the center of the grid)
	def add_detector(self, Npixel, lambda_min, lambda_max, N_lambda, alpha_x, alpha_y, D, dx=None, dy=None, delta_x=None, delta_y=None):
		self.Npixel_detector = Npixel
		if dx is None and dy is None and delta_x is None and delta_y is None:
				detector = (
					f'    <detector_dust nr_pixel = "{str(Npixel)}"> '
					f'{str(lambda_min)} {str(lambda_max)} {str(N_lambda)} 1 '
					f'{str(alpha_x)} {str(alpha_y)} {str(D)}'
					)		
				self.detectors.append(detector)
				return

		else:
			if delta_x is None and delta_y is None:
				detector = (
					f'    <detector_dust nr_pixel = "{str(Npixel)}"> '
					f'{str(lambda_min)} {str(lambda_max)} {str(N_lambda)} 1 '
					f'{str(alpha_x)} {str(alpha_y)} {str(D)} {str(dx)} {str(dy)}'
					)		
				self.detectors.append(detector)
				return
			else:
				detector = (
					f'    <detector_dust nr_pixel = "{str(Npixel)}"> '
					f'{str(lambda_min)} {str(lambda_max)} {str(N_lambda)} 1 '
					f'{str(alpha_x)} {str(alpha_y)} {str(D)} {str(dx)} {str(dy)} {str(delta_x)} {str(delta_y)}'
					)		
				self.detectors.append(detector)
				return
	
    # Add a spherical array of detectors 
	# nside: integer, and number of detectors = 12 * nside * nside
    # saved_path: the path to the /detector_angles/ folder which contains the POLARIS angles for the spherical array of detectors
    #             if left empty, it will calculate the angles using spherical_detectors (which is based on healpy)
	# Npixel, lambda_min, lambda_max, N_lambda, D, dx, dy, delta_x, delta_y: same as in add_detector
	def add_spherical_arrays(self, nside, Npixel, lambda_min, lambda_max, N_lambda, D, saved_path=None, dx=None, dy=None, delta_x=None, delta_y=None):
		self.detectors = []
		self.nside_spherical = nside
		if saved_path is None:
			theta, phi, alpha, beta = save_angle_detectors(saved_path, nside, save_coord=False, save_polaris_angle=False)
		else:
			alpha, beta = np.loadtxt(saved_path+f'/detector_angles/rotation_angles_detectors_nside{nside}.txt')

		for index_detector in range(len(alpha)):
			self.add_detector(Npixel=Npixel, lambda_min=lambda_min, lambda_max=lambda_max, 
							N_lambda=N_lambda, 
							alpha_x=np.round(alpha[index_detector]*180/np.pi, 3) % 360,
							alpha_y=np.round(beta[index_detector]*180/np.pi,3) % 360, 
							D=D,
							dx=dx, dy=dy,
							delta_x=delta_x, delta_y=delta_y
							)
		return
	
    
    

	# "path" is the path to a single dust parameters file, 
	# Ξi is the mass fraction of the material, 
	# q is the exponent of the grain size power-law distribution 
	# Nd(a) ∝ aq, and amin and amax are the dust grain radii wich 
	# have to be in the range as defined in the dust parameters file.
	def add_dust_component(self, path, xi, q, a_min, a_max):
		dust = "    <dust_component>     \""+path+"\" "+xi+" "+q+" "+a_min+" "+a_max
		self.dust_components.append(dust)
		return



	# Import massive stars from a torch snapshot and set them as stellar sources for POLARIS
	# particle_file: path to the torch particle file
	# minimum_mass: minimum mass of stars to be included (in solar masses)
	# num_photons: number of photons to be emitted by each star
	# 			   can be a single integer/float or a tuple for interpolation based on luminosity
	# stellar_metallicity: metallicity of the stars (default is 0.014)
	# spectra: type of spectra to use for the stars
	# 		   set to 'blackbody' for a blackbody spectrum
	# 		   set to 'ZAMS' for zero-age main sequence spectra from a precomputed SED grid
	# dir_sed_grid: directory containing the precomputed SED grid (only used if spectra is set to 'ZAMS')
	# dir_save_sed: directory to save the SED files for the stars (only used if spectra is set to 'ZAMS')
	def set_stellar_sources_from_torch(self, particle_file, minimum_mass, num_photons=1e6, stellar_metallicity=0.014, 
									   spectra='blackbody', dir_sed_grid='./', dir_save_sed=None):
		ds = yt.load(particle_file)
		ad = ds.all_data()

		# Find stars (which has particle_csgm == 0.0) that are more massive than the minimum mass
		star_idx = np.logical_and(ad['all', 'particle_csgm'] == 0.0, (ad['all', 'particle_old_pmass']*yt.units.g).to("Msun").v >= minimum_mass)
		stars = Particles(len(ad['all','particle_mass'][star_idx]))
		stars.initial_mass = (ad['all', 'particle_old_pmass'][star_idx]*yt.units.g).to("Msun").v | units.MSun
		stars.stellar_type = np.ones(len(stars)) | units.stellar_type

		# stellar positions in meters
		stars_x = ad['all', 'particle_posx'][star_idx].to('m').v
		stars_y = ad['all', 'particle_posy'][star_idx].to('m').v
		stars_z = ad['all', 'particle_posz'][star_idx].to('m').v
		
		# Calculate the stellar age in Myr based on the current simulation time and the particle creation time
		t_evol  = ds.current_time.in_units('Myr').v - (ad['all', 'particle_creation_time'][star_idx]*yt.units.s).to('Myr').v

		se = SeBa()
		se.initialize_code()

		# Evolve the star with SeBa based on stellar metallicity and age
		_tmp = se.evolve_star(stars.initial_mass, t_evol | units.Myr, stellar_metallicity) 
		se_time, se_mass, se_radius, se_lum, se_temp, se_evol_time, se_type = _tmp

		
		stars_r = se_radius.value_in(units.RSun) # in solar radii
		stars_T = se_temp.value_in(units.K) # in Kelvin
		stars_L = se_lum.value_in(units.LSun) # in solar luminosity
		stars_log_L = np.log10(stars_L)

		# Determine the number of photons for each star
		# If num_photons is a tuple, interpolate the number of photons based on the star's luminosity (in log space).
		if isinstance(num_photons, tuple):
			log_num_photons = np.interp(
				stars_log_L,
				(stars_log_L.min(), stars_log_L.max()),
				(np.log10(num_photons[0]), np.log10(num_photons[1])),
			)

			list_num_photons = np.rint(10**log_num_photons).astype(int)
		# If num_photons is a single integer or float, use that value for all stars.
		elif isinstance(num_photons, (int, float)):
			list_num_photons = np.full(len(stars_L), num_photons, dtype=int)
		else:
			raise ValueError("num_photons must be either a tuple or a single integer/float.")


		# Use blackbody spectra for the stars if spectra is set to 'blackbody', otherwise use ZAMS spectra from the SED grid
		if spectra == 'blackbody':
			for n in range(len(stars_x)):
				self.sources_string += "  <source_star nr_photons = \""+str(list_num_photons[n])+"\">	"+str(stars_x[n])+"	"+str(stars_y[n])+"	"+str(stars_z[n])+"	"+str(stars_r[n])+" "+str(stars_T[n])+"\n"
		
		# Use (scaled) zero-age main sequence (ZAMS) spectra for the stars if spectra is set to 'ZAMS'
		elif spectra == 'ZAMS':
			vec_M0 = stars.initial_mass.value_in(units.MSun)
			# Load the stellar SED grid from the specified directory
			stellar_sed_grid = np.load(dir_sed_grid+'stellar_sed_grid.npz')

			mass_Msun = stellar_sed_grid['mass_Msun']
			mets = stellar_sed_grid['mets']
			lams_A = stellar_sed_grid['lams_A']
			SEDs_Lsun_per_A_sorted = stellar_sed_grid['SEDs_Lsun_per_A_sorted']

			if stellar_metallicity < mets[0] or stellar_metallicity > mets[-1]:
				raise ValueError(f"Stellar metallicity {stellar_metallicity} is outside the range of the SED grid ({mets[0]} to {mets[-1]}). Please choose a metallicity within this range.")

			logM_grid = np.log10(mass_Msun)
			logZ_grid = np.log10(mets)

			log_SEDs_Lsun_per_A = np.log10(SEDs_Lsun_per_A_sorted)

			interp_log_SED = RegularGridInterpolator(
				(logM_grid, logZ_grid),
				log_SEDs_Lsun_per_A,
				method="linear",
				bounds_error=True
			)

			# Create a directory to save the SED files for the stars if it doesn't exist
			if dir_save_sed is not None:
				if not os.path.isdir(dir_save_sed+f'/SED_saved_Mmin{minimum_mass}'):
					os.mkdir(dir_save_sed+f'/SED_saved_Mmin{minimum_mass}')


			# Check if the star SED already exists by comparing the star properties with the saved star statistics
			write_SED_files = True
			stats_path = (
				dir_save_sed
				+ f"/SED_saved_Mmin{minimum_mass}/star_statistics.txt"
			)
			if os.path.isfile(stats_path):
				star_stats = np.loadtxt(stats_path, ndmin=2)
				current_star_stats = np.column_stack([
					stars_x,
					stars_y,
					stars_z,
					vec_M0,
					stars_L,
					list_num_photons,
				])
				# Check if the number of stars and their properties are the same as the saved statistics
				same_shape = star_stats.shape == current_star_stats.shape
				same_properties = (
					same_shape
					and np.allclose(
						star_stats,
						current_star_stats,
						rtol=1e-3,
						atol=0,
					)
				)
				if same_properties:
					write_SED_files = False

			# Write the SED files for each star if they don't already exist or if the properties have changed
			if write_SED_files:
				for index_star in range(len(stars_x)):
					# Interpolate the SED for the star's mass and metallicity (in log space) using the SED grid
					# Use maximum mass in the SED grid if the star's mass is above the maximum mass in the grid
					if vec_M0[index_star] > mass_Msun[-1]:
						log_SED_interpolated = interp_log_SED(np.array([[np.log10(mass_Msun[-1]), np.log10(stellar_metallicity)]]))[0]					
					else:
						log_SED_interpolated = interp_log_SED(np.array([[np.log10(vec_M0[index_star]), np.log10(stellar_metallicity)]]))[0]
					
					SED_interpolated = 10**log_SED_interpolated

					# Scale the SED to the luminosity computed by SeBa
					f_scale_SED = stars_L[index_star]/np.trapezoid(SED_interpolated, lams_A)

					if vec_M0[index_star] > mass_Msun[-1]:
						warnings.warn(f"Star {index_star} has mass {vec_M0[index_star]} Msun, which is above the maximum mass in the SED grid ({mass_Msun[-1]} Msun). Scaling the maximum mass SED by a factor of {f_scale_SED:.4f} for this star.")


					SED_interpolated = 10**log_SED_interpolated
					mask_nonzero = (log_SED_interpolated > -30)
					# POLARIS requires wavelength in meters and spectral luminosity in W/m, so we convert the units accordingly
					wavelength_save = lams_A[mask_nonzero]*1e-10 #1e-10 converts Angstroms to meters
					sed_save = f_scale_SED*SED_interpolated[mask_nonzero]*L_sun*1e10 #L_sun*1e10 converts Lsun/Angstrom to W/m
					data = np.column_stack([wavelength_save, sed_save])

					# Save the SED to a file for POLARIS input
					np.savetxt(
						dir_save_sed+f'/SED_saved_Mmin{minimum_mass}/SED_star{index_star}.dat',
						data,
						header="lambda [m]    L_lambda [W/m]",
						comments="# ",
						fmt="%.8e"
					)

					self.sources_string += "  <source_star nr_photons = \""+str(list_num_photons[index_star])+"\">	"+str(stars_x[index_star])+"	"+str(stars_y[index_star])+"	"+str(stars_z[index_star])+"	\""+dir_save_sed+f'/SED_saved_Mmin{minimum_mass}/SED_star{index_star}.dat'+"\"\n"
				
				np.savetxt(dir_save_sed+f'/SED_saved_Mmin{minimum_mass}/star_statistics.txt',
								np.column_stack([stars_x, stars_y, stars_z, vec_M0, stars_L, list_num_photons]),
								header="x [m]    y [m]    z [m]    M0 [Msun]    L [Lsun]    num_photons",)
			
			# If the SED files already exist and the properties haven't changed, just add the sources to the sources_string for POLARIS input
			else:
				print("SED files already exist and star properties haven't changed. Using existing SED files.")
				for index_star in range(len(stars_x)):
					if vec_M0[index_star] > mass_Msun[-1]:
						warnings.warn(f"Star {index_star} has mass {vec_M0[index_star]} Msun, which is above the maximum mass in the SED grid ({mass_Msun[-1]} Msun). Scaling the maximum mass SED for this star.")

					self.sources_string += "  <source_star nr_photons = \""+str(list_num_photons[index_star])+"\">	"+str(stars_x[index_star])+"	"+str(stars_y[index_star])+"	"+str(stars_z[index_star])+"	\""+dir_save_sed+f'/SED_saved_Mmin{minimum_mass}/SED_star{index_star}.dat'+"\"\n"

		return


	# Set stellar sources from list of star positions, radii, and temperatures (assuming blackbody spectra)
	# This function is deprecated
	def set_stellar_sources(self, x, y, z, rad, temp, num_photons="1e6"):
		for n in range(len(x)):
			self.sources_string += "  <source_star nr_photons = \""+num_photons+"\">	"+str(x[n])+"	"+str(y[n])+"	"+str(z[n])+"	"+str(rad[n])+" "+str(temp[n])+"\n"
		return



	# Generate the command file for POLARIS
	# task_name: "temp" for dust temperature calculation
	# 			 "em_stellar" for stellar+dust emission
	# 			 "em_dust" for dust only emission
	# command_filename: name of the command file to be generated
	# write_mode: 'a' for append, 'w' for write
	# sed_label: optional label for the SED output directory, renames the emission output directory to avoid overwriting previous results
	def generate_command_file(self, task_name, command_filename, write_mode='a', sed_label=None, dust_offset=False):
		cmd_file = open(command_filename, write_mode)

		if task_name == "temp":
			if not os.path.isdir(self.directory+"/temp"):
				os.mkdir(self.directory+"/temp")
			# CMD_TEMP runs a simulation for heating the dust by considering different photon emitting sources
			cmd_file.write("<task> 1\n  <cmd>  CMD_TEMP\n\n")
			for d_c in self.dust_components:
				cmd_file.write(d_c+"\n")

			# path to input grid
			cmd_file.write("  <path_grid>      \""+self.gridfile+"\"\n")
			# path for the temperature output files
			cmd_file.write("  <path_out>      \""+self.directory+"/temp/"+"\"\n")
			cmd_file.write("  <write_inp_midplanes> "+self.n_bins+"\n")
			cmd_file.write("  <write_out_midplanes> "+self.n_bins+"\n")
			cmd_file.write("  <conv_dens> "+self.conv_dens+"\n")
			cmd_file.write("  <conv_len> "+self.conv_len+"\n")
			cmd_file.write("  <mu> "+self.mu+"\n")
			if not (self.mass_fraction is None):
				cmd_file.write("  <mass_fraction> "+self.mass_fraction+"\n")
			if dust_offset:
				cmd_file.write("  <dust_offset> 1 \n")
			cmd_file.write("  <sub_dust> 1 \n")
			cmd_file.write("  <nr_threads> "+self.num_threads+"\n")
			cmd_file.write(self.sources_string+"\n")
			cmd_file.write("</task>"+"\n\n\n")

		elif task_name == "em_stellar":
			if self.nside_spherical is not None:
				if sed_label is None:
					dir_sed = self.directory+f"/em_stellar_dust_nside{int(self.nside_spherical)}"
				else:
					dir_sed = self.directory+f"/em_stellar_dust_nside{int(self.nside_spherical)}_{sed_label}"
			else:
				if sed_label is None:
					dir_sed = self.directory+f"/em_stellar_dust_Npixel{int(self.Npixel_detector)}"
				else:
					dir_sed = self.directory+f"/em_stellar_dust_Npixel{int(self.Npixel_detector)}_{sed_label}"

			if not os.path.isdir(dir_sed):
				os.mkdir(dir_sed)


			# CMD_DUST_EMISSION defines the ray-tracing
			# It can be run only after CMD_TEMP has been run, as it requires the dust temperature grid as input
			cmd_file.write("<task> 1\n  <cmd>  CMD_DUST_EMISSION\n\n")
			for d_c in self.dust_components:
				cmd_file.write(d_c+"\n")
			cmd_file.write("    <path_grid>      \""+self.directory+"/temp/grid_temp.dat"+"\"\n")
			# path for the  output files
			cmd_file.write("    <path_out>       \""+dir_sed+"/"+"\"\n")
			# The conversion factor is set to 1.0 because we are reading the grid from CMD_TEMP, which is already converted into the SI unit
			cmd_file.write("  <conv_dens>     1.0\n")  
			cmd_file.write("  <conv_len>      1.0\n") 
			cmd_file.write("  <mu> "+self.mu+"\n")

			# The dust-to-gas ratio, if the grid does not have dust densities. Mostly deprecated, as the grid should have dust densities.
			if not (self.mass_fraction is None):
				cmd_file.write("  <mass_fraction> "+self.mass_fraction+"\n")
			
			cmd_file.write("  <nr_threads> "+self.num_threads+"\n")
			cmd_file.write("  <max_subpixel_lvl> "+self.subpixel_lvl+"\n")
			# detectors
			for det in self.detectors:
				cmd_file.write(det+"\n")
			cmd_file.write(self.sources_string+"\n")
			cmd_file.write("</task>"+"\n\n\n")

		elif task_name == "em_dust":
			if self.nside_spherical is not None:
				dir_sed_dust = self.directory+f"/em_dust_nside{int(self.nside_spherical)}"
				if not os.path.isdir(dir_sed_dust):
					os.mkdir(self.directory+f"/em_dust_nside{int(self.nside_spherical)}")
			else:
				dir_sed_dust = self.directory+f"/em_dust_Npixel{int(self.Npixel_detector)}"
				if not os.path.isdir(dir_sed_dust):
					os.mkdir(self.directory+f"/em_dust_Npixel{int(self.Npixel_detector)}")

			# CMD_DUST_EMISSION defines the ray-tracing
			# It can be run only after CMD_TEMP has been run, as it requires the dust temperature grid as input
			cmd_file.write("<task> 1\n  <cmd>  CMD_DUST_EMISSION\n\n")
			for d_c in self.dust_components:
				cmd_file.write(d_c+"\n")
			cmd_file.write("    <path_grid>      \""+self.directory+"/temp/grid_temp.dat"+"\"\n")
			# path for the  output files
			cmd_file.write("  <path_out>      \""+dir_sed_dust+"/"+"\"\n")
			cmd_file.write("  <conv_dens>     1.0\n")
			cmd_file.write("  <conv_len>      1.0\n")
			cmd_file.write("  <mu> "+self.mu+"\n")
			if not (self.mass_fraction is None):
				cmd_file.write("  <mass_fraction> "+self.mass_fraction+"\n")
			cmd_file.write("  <nr_threads> "+self.num_threads+"\n")
			# detectors
			for det in self.detectors:
				cmd_file.write(det+"\n")
			cmd_file.write("</task>"+"\n\n\n")

			cmd_file.close()

		return


