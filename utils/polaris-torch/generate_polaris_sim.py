
# import function to convert FLASH grid
import flash_to_polaris
from Polaris_torch import Polaris

# To use full stellar spectra, set IGNORE_WAVELENGTH_RANGE true in 
# POLARIS/src/Typedefs.hpp and recompile with ./compile.sh -u 

# -------------------- INPUTS -------------------- #

torch_gridfile = "turbsph_hdf5_plt_cnt_0389"
torch_partfile = "turbsph_hdf5_part_0389"
sim_dir="./output"

# -------------------- CONVERT GRID -------------------- #

polaris_gridfile = torch_gridfile+"_polaris.dat"
flash_to_polaris.convert(torch_gridfile, polaris_gridfile)

# -------------------- CREATE POLARIS OBJECT -------------------- #
polaris = Polaris(sim_directory=sim_dir, gridfile=polaris_gridfile,
		  n_bins="256",conv_dens="1e6",conv_len="0.01", mu="1.3", 
                  mass_fraction="0.01", num_threads="32")

# -------------------- ADD DUST COMPONENTS -------------------- #
# "path" is the path to a single dust parameters file, 
# xi is the mass fraction of the material, 
# q is the exponent of the grain size power-law distribution 
# Nd(a) ∝ aq, and amin and amax are the dust grain radii wich 
# have to be in the range as defined in the dust parameters file.

polaris.add_dust_component(path="/home/bpolak/POLARIS/input/dust_cs/silicate_oblate.dat",
							xi="0.625", q="-3.5", a_min="5.0e-9", 
							a_max="250.0e-9")
polaris.add_dust_component(path="/home/bpolak/POLARIS/input/dust_cs/graphite_oblate.dat",
							xi="0.375", q="-3.5", a_min="5.0e-9", 
							a_max="250.0e-9")

# -------------------- ADD DETECTORS -------------------- #
# Npixel is the number of pixel, λmin is the shortest wavelength, 
# λmax is the longest wavelength,Nλ is the amount of wavelengths 
# used for simulation (if 1, only λmin is used), α_x and α_y are 
# the rotation angles around the first and second rotation axis, 
# and D is the distance to the observer.
# Default polaris wavelength range: 2.63e-07 [m] to 0.003 [m]
# See note at top of file to expand spectrum


# Define 6 plane detectors for each cube surface to obtain escape fractions
alpha_xs = [0,0,0,90,180,270]
alpha_ys = [0,90,-90,0,0,0]

for a in range(6):
    # Full spectrum
    polaris.add_detector(Npixel="256", lambda_min="2.63e-7", lambda_max="3e-3", 
     					 N_lambda="10", alpha_x=alpha_xs[a], alpha_y=alpha_ys, D="3.086e+18")

# -------------------- ADD SOURCES -------------------- #
# set minimum mass for tracing stellar radiation
polaris.set_stellar_sources_from_torch(torch_partfile, minimum_mass=20.0, num_photons="1e6")

# Can also define sources with an array of x,y,z,radius,temp of stars
# polaris.set_stellar_sources(x, y, z, rad, temp, num_photons="1e6")

# -------------------- GENERATE COMMAND FILE -------------------- #
polaris.generate_command_file(command_filename="torch_cmd")

# NOW you can just run ./polaris torch.cmd
