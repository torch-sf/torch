# ------- GENERATE POLARIS SIMULATION OF A TORCH SNAPSHOT ------- #
# Brooke Polak, January 2025
import os
import numpy as np
import yt
from amuse.lab import *
from amuse.units import units

# List of functions:
	# set_stellar_sources(self, particle_file)
	# set_stellar_sources(self, x, y, z, rad, temp)

# PARAMETERS:
# grid file
# output directory
# number of bins
# conv_dens 
# conv_len
# mu
# mass_fraction
# nr_threads
# list of stellar sources (x,y,z, radius, temperature)
# detectors X
# dust components X

class Polaris:

	# Class constructor
	def __init__(self, sim_directory, gridfile, n_bins="256", conv_dens="1e6",
				 conv_len="0.01", mu="1.3", mass_fraction="0.01", num_threads="-1", 
				 dust_components=[], detectors=[]):

		self.directory = sim_directory
		# make output directories
		if not os.path.isdir(self.directory):
			os.mkdir(self.directory)

		if not os.path.isdir(self.directory+"/temp"):
			os.mkdir(self.directory+"/temp")
			os.mkdir(self.directory+"/em_stellar_dust")
			os.mkdir(self.directory+"/em_dust")

		self.gridfile = gridfile
		self.n_bins = n_bins
		self.conv_dens = conv_dens
		self.conv_len = conv_len
		self.mu = mu
		self.mass_fraction = mass_fraction
		self.num_threads = num_threads
		self.dust_components = []
		self.detectors = []
		# create chunk of sources string once to avoid multiple for loops
		self.sources_string = ""

		return 

	# Npixel is the number of pixel, λmin is the shortest wavelength, 
	# λmax is the longest wavelength,Nλ is the amount of wavelengths 
	# used for simulation (if 1, only λmin is used), α_x and α_y are 
	# the rotation angles around the first and second rotation axis, 
	# and D is the distance to the observer.
	def add_detector(self, Npixel, lambda_min, lambda_max, N_lambda, alpha_x, alpha_y, D):
		detector = "    <detector_dust nr_pixel = \""+Npixel+"\"> "+lambda_min+" "+lambda_max+" "+N_lambda+" 1 "+alpha_x+" "+alpha_y+" "+D
		self.detectors.append(detector)
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

	# Use SeBa to get the temperature and radii of evolved stars
	def set_stellar_sources_from_torch(self, particle_file, minimum_mass, num_photons="1e6"):
		ds = yt.load(particle_file)
		ad = ds.all_data()

		### Get the total UV energy emitted for a given star mass
		se = SeBa()
		se.initialize_code()

		star_idx = np.logical_and(ad['all', 'particle_csgm'] == 0.0, (ad['all', 'particle_old_pmass']*yt.units.g).to("Msun").v >= minimum_mass)
		stars = Particles(len(ad['all','particle_mass'][star_idx]))
		stars.initial_mass = (ad['all', 'particle_old_pmass'][star_idx]*yt.units.g).to("Msun").v | units.MSun
		stars.stellar_type = np.ones(len(stars)) | units.stellar_type
		t_evol  = ds.current_time.in_units('Myr').v - (ad['all', 'particle_creation_time'][star_idx]*yt.units.s).to('Myr').v

		# metallicity 0.02
		_tmp = se.evolve_star(stars.initial_mass, t_evol | units.Myr, 0.02) 
		se_time, se_mass, se_radius, se_lum, se_temp, se_evol_time, se_type = _tmp

		# cm 
		stars_x = ad['all', 'particle_posx'][star_idx].to('m').v
		stars_y = ad['all', 'particle_posy'][star_idx].to('m').v
		stars_z = ad['all', 'particle_posz'][star_idx].to('m').v
		# solar radii
		stars_r = se_radius.value_in(units.RSun)
		# Kelvin
		stars_T = se_temp.value_in(units.K)

		for n in range(len(stars_x)):
			self.sources_string += "  <source_star nr_photons = \""+num_photons+"\">	"+str(stars_x[n])+"	"+str(stars_y[n])+"	"+str(stars_z[n])+"	"+str(stars_r[n])+" "+str(stars_T[n])+"\n"

		return

	# Set stellar sources from list of star positions, radii, and temperatures
	def set_stellar_sources(self, x, y, z, rad, temp, num_photons="1e6"):
		for n in range(len(x)):
			self.sources_string += "  <source_star nr_photons = \""+num_photons+"\">	"+str(x[n])+"	"+str(y[n])+"	"+str(z[n])+"	"+str(rad[n])+" "+str(temp[n])+"\n"
		return

	def generate_command_file(self, command_filename):

		cmd_file = open(command_filename, 'w')
		# CMD_TEMP runs a simulation for heating the dust by considering 
		# different photon emitting sources
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
		cmd_file.write("  <mass_fraction> "+self.mass_fraction+"\n")
		cmd_file.write("  <num_threads> "+self.num_threads+"\n")
		cmd_file.write(self.sources_string+"\n")
		cmd_file.write("</task>"+"\n\n\n")

		# CMD_DUST_EMISSION defines the ray-tracing
		cmd_file.write("<task> 1\n  <cmd>  CMD_DUST_EMISSION\n\n")
		for d_c in self.dust_components:
			cmd_file.write(d_c+"\n")
		cmd_file.write("    <path_grid>      \""+self.directory+"/temp/grid_temp.dat"+"\"\n")
		# path for the  output files
		cmd_file.write("    <path_out>       \""+self.directory+"/em_stellar_dust/"+"\"\n")
		cmd_file.write("  <conv_dens>     1.0\n")
		cmd_file.write("  <conv_len>      1.0\n")
		cmd_file.write("  <mu> "+self.mu+"\n")
		cmd_file.write("  <mass_fraction> "+self.mass_fraction+"\n")
		cmd_file.write("  <num_threads> "+self.num_threads+"\n")
		# detectors
		for det in self.detectors:
			cmd_file.write(det+"\n")
		cmd_file.write(self.sources_string+"\n")
		cmd_file.write("</task>"+"\n\n\n")

		# CMD_DUST_EMISSION defines the ray-tracing
		cmd_file.write("<task> 1\n  <cmd>  CMD_DUST_EMISSION\n\n")
		for d_c in self.dust_components:
			cmd_file.write(d_c+"\n")
		cmd_file.write("    <path_grid>      \""+self.directory+"/temp/grid_temp.dat"+"\"\n")
		# path for the  output files
		cmd_file.write("  <path_out>      \""+self.directory+"/em_dust/"+"\"\n")
		cmd_file.write("  <conv_dens>     1.0\n")
		cmd_file.write("  <conv_len>      1.0\n")
		cmd_file.write("  <mu> "+self.mu+"\n")
		cmd_file.write("  <mass_fraction> "+self.mass_fraction+"\n")
		cmd_file.write("  <num_threads> "+self.num_threads+"\n")
		# detectors
		for det in self.detectors:
			cmd_file.write(det+"\n")
		cmd_file.write("</task>"+"\n\n\n")

		cmd_file.close()

		return




