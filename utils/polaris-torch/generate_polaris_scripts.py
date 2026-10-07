# First Written by Brooke Polak, January 2025
# Edited by Donglin Wu, July 2026


import flash_to_polaris
from Polaris_torch import Polaris
import subprocess
import numpy as np
import sys

print('Modules loaded')

# To use full stellar spectra, set IGNORE_WAVELENGTH_RANGE true in 
# POLARIS/src/Typedefs.hpp and recompile with ./compile.sh -u 



# Main function to generate POLARIS command file from torch simulations
#
# -------------------- POLARIS setup -------------------- #
# run_mode: "T_dust", create a POLARIS grid from Torch files and run POLARIS to compute the dust temperature
#           "detectors_spherical", run POLARIS with existing grid to derive emission from spherical detectors
#                                  primary purpose is to obtain FUV escape fraction
#           "detectors_image", run POLARIS existing grid to obtain simulated observations from six detectors (six faces of a cube)
#                              images are the desired outputs
# dir_polaris: path to the POLARIS installation directory
#
# -------------------- PHYSICAL STATE OF THE SYSTEM -------------------- #
# r_dust_to_gas: dust-to-gas ratio
# stellar_metallicity: metallicity of the stellar population
# minimum_mass_stars: mass of stars considered 
# mu: mean molecular weight, the default is 1.3
#     we generally ignore any dependence of mean molecular weight on metallicity
#
# -------------------- TORCH & OTHER INPUT -------------------- #
# snapshot_number: the index of the snapshot
# dir_torch: path to the directory containing the Torch output files
# dir_sed_grid: path to the directory containing the SED grid files for stellar sources
# 
# -------------------- POLARIS OUTPUT SETUP -------------------- #
# dir_output: path to the directory where POLARIS output will be saved
# extra_label: additional label to distinguish different runs with different grid
# sed_label: additional label to distinguish different SED runs using the same grid
#
# -------------------- POLARIS OUTPUT DETECTORS -------------------- #
# Two options for detector wavelengths: 1) input a list of observed wavelengths (in meters) and redshift, 
#                                       2) input a wavelength range and number of wavelengths
# 1) a list of observed wavelengths (in meters)
# wavelength_bands: observed wavelengths for the detectors on Earth
# redshift: redshift of the system, used to convert observed wavelengths to rest-frame wavelengths
# 2) a wavelength range and number of wavelengths
# wl_min: minimum wavelength (in meters) for the detectors
# wl_max: maximum wavelength (in meters) for the detectors
# N_lambda: number of wavelengths for the detectors
# 
# d_to_source: distance from the detector to the center of the grid in meters
#
# if run_mode is set to "detectors_spherical":
#       nside_spherical: nside for healpy to generate detectors, number of detectors = 12*nside*nside
#       dir_detector_angles_saved: the path to the /detector_angles/ folder which contains the POLARIS angles for the spherical array of detectors
# if run_mode is set to "detectors_image":
#       N_pixel_image: number of pixels in each detector image, default is 256
#
#
# -------------------- OTHER PARAMETERS THAT ARE NOT ARGUMENTS -------------------- #
# dust components: the dust components used in POLARIS, currently set to the default THEMIS model
# dust grain size distribution: currently set to a power-law distribution with exponent -3.5 and bounds 5 nm to 250 nm
# threshold for sputtering timescale: default is to use the current simulation time as the threshold for sputtering timescale 
# level of subpixeling used by POLARIS: default is 3, recommended to do a convergence test on this parameter for infrared




def generate_polaris_cmd(run_mode, r_dust_to_gas, stellar_metallicity, 
                         snapshot_number, 
                         dir_polaris, dir_torch, dir_output,
                         extra_label,
                         minimum_mass_stars=20.0,
                         mu=1.3,
                         dir_sed_grid=None, 
                         sed_label=None,
                         nside_spherical=None,
                         dir_detector_angles_saved=None,
                         N_pixel_image=256,
                         d_to_source=3.086e18,
                         wavelength_bands=None, redshift=None, 
                         wl_min=None, wl_max=None, N_lambda=None):

    label_snapshot = f"{int(snapshot_number):04d}"

    # Torch particle and grid files
    # Change the following if the Torch output file naming convention is different
    torch_gridfile = f"/turbsph_hdf5_plt_cnt_{label_snapshot}"
    torch_partfile = f"/turbsph_hdf5_part_{label_snapshot}"

    if N_lambda is None:
        N_lambda = len(wavelength_bands)
    sim_dir=dir_output+f'/output_{label_snapshot}_Nlambda{N_lambda}_{extra_label}'


    # -------------------- CONVERT GRID -------------------- #
    polaris_gridfile = dir_output+torch_gridfile+f"_polaris_{label_snapshot}_{extra_label}.dat"
    if run_mode == "T_dust":
        flash_to_polaris.convert(dir_torch+torch_gridfile, polaris_gridfile, 
                                 mass_fraction=r_dust_to_gas, 
                                 mu=mu,
                                 # Default: use the current simulation time as the threshold for sputtering timescale.
                                 # If you want to use a fixed threshold, set this to a float value in Gyr.
                                 # If you want to disable sputtering, set this to a very small value (e.g., 1e-30). 
                                 timescale_sputtering="time_snapshot",
                                 # Set to True to clear dust in cells within the sublimation radius of any star.
                                 # For Torch runs, this is unnecessary, as the grid cell sizes are much larger than the sublimation radius of any star.
                                 clear_sublimation_regions=False, 
                                 file_particle_input=dir_torch+torch_partfile,
                                 minimum_mass_stars=minimum_mass_stars,
                                 stellar_metallicity=stellar_metallicity,
                                 # Sublimation temperature of the dust
                                 # Should be consistent with the dust components used (see below)
                                 T_sub=2000.0)


    # -------------------- CREATE POLARIS OBJECT -------------------- #
    polaris = Polaris(sim_directory=sim_dir, gridfile=polaris_gridfile,
            n_bins="256",conv_dens="1e6",conv_len="0.01", 
            mu=f"{mu}", 
            mass_fraction=None,
            # Default: use all processors available 
            num_threads="-1", 
            # Level of subpixeling used by POLARIS
            # For infrared, it is recommended to do a convergence test on this parameter.
            subpixel_lvl="3")

    # -------------------- ADD DUST COMPONENTS -------------------- #
    # path: the path to the dust parameters files
    #       if using the default THEMIS model, they should be part of the POLARIS directory as specified below 
    # xi: the mass fraction of the material in the dust mixture, where the sum of all xi should be 1.0
    # q, amin and amax: the exponent and bounds of the grain size power-law distribution 
    # WARNING: if the grain size distribution is changed, 
    #          the mass-weighted average of the grain size should be changed in dust_helper.py

    polaris.add_dust_component(path=dir_polaris+"/input/dust_cs/aOlM5.dat",
                                xi="0.45", q="-3.5", a_min="5.0e-9", 
                                a_max="250.0e-9")
    polaris.add_dust_component(path=dir_polaris+"/input/dust_cs/aPyM5.dat",
                                xi="0.15", q="-3.5", a_min="5.0e-9", 
                                a_max="250.0e-9")
    polaris.add_dust_component(path=dir_polaris+"/input/dust_cs/CM20.dat",
                                xi="0.4", q="-3.5", a_min="5.0e-9",
                                a_max="250.0e-9")

    # -------------------- ADD SOURCES -------------------- #
    # Include sources from Torch particle data.
	# num_photons: number of photons to be emitted by each star
	# 			   can be a single integer/float or a tuple for interpolation based on luminosity
    polaris.set_stellar_sources_from_torch(dir_torch+torch_partfile, minimum_mass=minimum_mass_stars, num_photons=(1e4, 1e6),
                                           stellar_metallicity=stellar_metallicity, spectra='ZAMS', dir_sed_grid=dir_sed_grid,
                                           dir_save_sed=sim_dir)




    if run_mode == "T_dust":
        polaris.generate_command_file("temp", write_mode='w', command_filename=dir_output+f'/torch_cmd_uv_{label_snapshot}_{extra_label}')


    # -------------------- ADD DETECTORS -------------------- #
    if run_mode == "detectors_spherical":
        if wavelength_bands is not None:
            raise ValueError("For spherical detector geometry, please specify a wavelength range (wl_min, wl_max) and number of wavelengths (N_lambda), instead of a list of wavelengths.")
        if nside_spherical is None:
            raise ValueError("For spherical detector geometry, please specify nside_spherical.")
            
        polaris.add_spherical_arrays(
                nside=nside_spherical,
                # Number of pixels in each detector image.
                # For UV, the default is 1 pixel, as we are only interested in the total flux.
                # For IR, a convergence test on this parameter (complementary to level of subpixeling) is recommended.
                Npixel="1", 
                lambda_min=f"{wl_min:.4e}", lambda_max=f"{wl_max:.4e}", N_lambda=str(N_lambda),
                # Distance to the source in m
                # Set to 100 pc for spherical detectors by default
                D=f"{d_to_source:.4e}",
                saved_path=dir_detector_angles_saved,
                )
        
        if sed_label is None:
            polaris.generate_command_file("em_stellar", write_mode='w', command_filename=dir_output+f'/torch_cmd_uv_{label_snapshot}_{extra_label}_spherical')
        else:
            polaris.generate_command_file("em_stellar", write_mode='w', command_filename=dir_output+f'/torch_cmd_uv_{label_snapshot}_{extra_label}_spherical_{sed_label}', sed_label=sed_label)



    elif run_mode == "detectors_image":

        # Angles at which the detectors are placed, in degrees.
        # Default for taking images is on six faces of a cube, with the cloud at the center.
        alpha_xs = [0,0,0,90,180,270]
        alpha_ys = [0,90,-90,0,0,0]

        if wavelength_bands is not None:
            if redshift is not None:
                wavelength_restframe = wavelength_bands/(1+redshift)
            else:
                wavelength_restframe = wavelength_bands

            for a in range(len(alpha_xs)):
                for index_wl, wl in enumerate(wavelength_restframe):
                    polaris.add_detector(Npixel=f"{N_pixel_image}", 
                                        lambda_min=f"{wl:.4e}", lambda_max=f"{wl:.4e}", N_lambda="1",
                                        alpha_x=alpha_xs[a], alpha_y=alpha_ys[a], D=f"{d_to_source:.4e}",)

        else:
            for a in range(len(alpha_xs)):
                polaris.add_detector(Npixel=f"{N_pixel_image}", 
                                    lambda_min=f"{wl_min:.4e}", lambda_max=f"{wl_max:.4e}", N_lambda=str(N_lambda),
                                    alpha_x=alpha_xs[a], alpha_y=alpha_ys[a], D=f"{d_to_source:.4e}",)


        if sed_label is None:
            polaris.generate_command_file("em_stellar", write_mode='w', command_filename=dir_output+f'/torch_cmd_uv_{label_snapshot}_{extra_label}_images')
        else:
            polaris.generate_command_file("em_stellar", write_mode='w', command_filename=dir_output+f'/torch_cmd_uv_{label_snapshot}_{extra_label}_images_{sed_label}', sed_label=sed_label)





