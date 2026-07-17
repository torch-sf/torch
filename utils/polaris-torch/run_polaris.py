import subprocess
import traceback
import numpy as np
import sys
import argparse

import os
from generate_polaris_scripts import generate_polaris_cmd

# Using subprocess to run commands in terminal
def run_command(cmd, workdir, timeout=None):
    result = subprocess.run(
        cmd,
        cwd=workdir,
        text=True,
        shell=True,
        timeout=timeout,
        check=False,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"Command {cmd!r} failed with code {result.returncode}"
        )


# Main function to generate POLARIS script and run POLARIS
# For the parameters, please refer to generate_polaris_scripts.py for more details.
def run_polaris(run_mode, r_dust_to_gas, stellar_metallicity, 
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

    generate_polaris_cmd(run_mode, r_dust_to_gas, stellar_metallicity, 
                         snapshot_number, 
                         dir_polaris, dir_torch, dir_output,
                         extra_label,
                         minimum_mass_stars=minimum_mass_stars,
                         mu=mu,
                         dir_sed_grid=dir_sed_grid, 
                         sed_label=sed_label,
                         nside_spherical=nside_spherical,
                         dir_detector_angles_saved=dir_detector_angles_saved,
                         N_pixel_image=N_pixel_image,
                         d_to_source=d_to_source,
                         wavelength_bands=wavelength_bands, redshift=redshift, 
                         wl_min=wl_min, wl_max=wl_max, N_lambda=N_lambda)
    
    polaris_exec = "polaris"

    if run_mode == "T_dust":
        print("Running POLARIS for dust temperature calculation...")
        run_command(f"{polaris_exec} {dir_output}/torch_cmd_uv_{label_snapshot}_{extra_label}", dir_polaris)
    
    elif run_mode == "detectors_images":
        print("Running POLARIS for detector images...")
        if sed_label is None:
            run_command(f"{polaris_exec} {dir_output}/torch_cmd_uv_{label_snapshot}_{extra_label}_images", dir_polaris)
        else:
            run_command(f"{polaris_exec} {dir_output}/torch_cmd_uv_{label_snapshot}_{extra_label}_images_{sed_label}", dir_polaris)
    
    elif run_mode == "detectors_spherical":
        print("Running POLARIS for an array of spherical detectors...")
        if sed_label is None:
            run_command(f"{polaris_exec} {dir_output}/torch_cmd_uv_{label_snapshot}_{extra_label}_spherical", dir_polaris)
        else:
            run_command(f"{polaris_exec} {dir_output}/torch_cmd_uv_{label_snapshot}_{extra_label}_spherical_{sed_label}", dir_polaris)



parser = argparse.ArgumentParser()
parser.add_argument("-r", "--run_mode", default='T_dust', required=True, type=str, 
                    help="The mode in which to run POLARIS (e.g., 'T_dust', 'detectors_images', 'detectors_spherical').")
parser.add_argument("-p", "--path_snapshot", default='./', required=False, type=str, 
                    help="Path to the Torch snapshot file.")
parser.add_argument("-f", "--file_snapshot", default='', required=False, type=str, 
                    help="Path to the file containing the list of snapshots to run POLARIS on. Each line should contain a comma-separated list of snapshot numbers.")
parser.add_argument("-l", "--line_snapshot", default=0, required=False, type=int, 
                    help="The line number in the snapshot list file for POLARIS to process.")
parser.add_argument("-s", "--snapshot_number", default=None, required=False, type=int, 
                    help="The snapshot number to run POLARIS on.")
parser.add_argument("-o", "--output_directory", default='', required=True, type=str, 
                    help="The directory where POLARIS output will be saved.")
parser.add_argument("-ep", "--epsilon", default=1e-6, required=False, type=float, 
                    help="The dust to gas ratio to use in POLARIS.")
parser.add_argument("-ns", "--nside", default=7, required=False, type=int, 
                    help="The nside for the spherical array of detectors to use in POLARIS.")
parser.add_argument("-el", "--extra_label", default='run1', required=False, type=str, 
                    help="The extra label to use in the POLARIS command file name.")
parser.add_argument("-sl", "--sed_label", default=None, required=False, type=str, 
                    help="The sed label to use in the POLARIS command file name.")
args = parser.parse_args()

# Determine which run mode to use for POLARIS
run_mode = args.run_mode

# Determine the list of snapshots to run POLARIS on
if args.snapshot_number is None:
    try:
        array_snapshot = []
        with open(args.file_snapshot, "r") as f:
            for line in f:
                arr = [int(x) for x in line.split(',')]
                array_snapshot.append(arr)    
        list_snapshot = array_snapshot[args.line_snapshot]
    except:
        raise ValueError("Please provide a valid file containing the list of snapshots and a valid line number.")
else:
    list_snapshot = [args.snapshot_number]
print('Snapshots running:', list_snapshot)


# Determine the output directory and create it if it doesn't exist
dir_output = args.output_directory
if not os.path.exists(dir_output):
    os.makedirs(dir_output)

# Determine the extra label and sed label for the POLARIS command file name
extra_label = args.extra_label
sed_label = args.sed_label

# Determine the dust to gas ratio and nside for spherical detectors
r_dust_to_gas = args.epsilon
nside_spherical = args.nside
if run_mode == "detectors_spherical":
    if nside_spherical is None or nside_spherical <= 0 or nside_spherical > 20:
        raise ValueError("Please provide a valid nside for the spherical array of detectors.")



# Define the directories for POLARIS, Torch output, and the current run
dir_polaris = "/home/dwu/mendel-nas1/POLARIS/POLARIS/"
dir_torch = args.path_snapshot
dir_run = "./"



for snapshot_number in list_snapshot:
    print(f'Running snapshot {snapshot_number}')
    try:
        run_polaris(run_mode, r_dust_to_gas, 0.01*0.014, 
                snapshot_number, 
                dir_polaris, dir_torch, dir_output,
                extra_label,
                minimum_mass_stars=20.0,
                mu=1.3,
                dir_sed_grid=dir_run, 
                sed_label=sed_label,
                nside_spherical=nside_spherical,
                dir_detector_angles_saved=None,
                N_pixel_image=256,
                d_to_source=3.086e18,
                wl_min=0.092*1e-6, wl_max=0.7*1e-6, N_lambda=40)
        print(f'Snapshot {snapshot_number} finished.')
    except Exception as e:
        print(e)
        traceback.print_exc()
        print(f'Snapshot {snapshot_number} failed.')


