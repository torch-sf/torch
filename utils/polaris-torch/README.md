## Utility to run Polaris on a Torch simulation snapshot
First written by Stefan Reissl (flash-to-polaris.py) and Brooke Polak (Polaris-torch.py, generate\_polaris\_scripts.py)

Edited by Donglin Wu


### Step 1: Download and install Polaris
Download and install Polaris from https://github.com/polaris-MCRT/POLARIS


### Step 2: Install dependencies
This tool does not require a full installation of Torch. However, the following dependencies are required:
- **Common dependencies**: `numpy`, `scipy`, `argparse` (`argparse` is included in the Python standard library)
- **Dependencies inherited from Torch**: `yt`, `amuse`

    Note that only `amuse-seba` and `amuse-framework` are required for this tool. Installation instructions for `amuse` are available at: https://amuse.readthedocs.io/en/latest/install/installing.html
- **Dependencies specific to this tool**: `healpy`

### Step 3: Customize the Python scripts for your application
Edit script generate\_polaris\_scripts.py. This python file includes a function to produce command files for POLARIS. 

* **IMPORTANT**: The default setting is to use all processors available, which works well for a slurm job that specifies the number of processors to be used. However, if the tool is used without a slurm job, please change line 128 of generate\_polaris\_scripts.py and specify the number of processors to be used, e.g., **num_threads="16"**.
* This tool assumes the default Torch output naming convention: "turbsph_hdf5_plt_cnt_0001" for grid files and "turbsph_hdf5_part_0001" for particle files. If a different naming convention is used, lines 93 and 94 of generate\_polaris\_scripts.py have to be modified. 
* The default dust model used is the [THEMIS model](https://www.ias.u-psud.fr/themis/THEMIS_model.html), which is built-in in the standard POLARIS installation. If different dust components should be used, relevant sections (starting line 133) in generate\_polaris\_scripts.py have to be modified. 
* Besides the changes described above, other modifications are more specialized and are documented in the comments within generate\_polaris\_scripts.py.

Edit script run\_polaris.py. This python file uses generate\_polaris\_scripts.py to write command files for POLARIS and runs terminal command of POLARIS using subprocess. 

* This script (starting line 152) is primarily where parameters can be adjusted. The parameters are summarized in the following table:

    | Parameter | Location | Default |
    |---|---|---|
    | Stellar metallicity | Line 154 | 1.4e-4 (0.01 solar) |
    | Minimum mass of stars to be included in POLARIS in solar mass | Line 160 | 20 |
    | Mean molecular weight | Line 161 | 1.3 |
    | Distance between the detector and center of the grid in meters | Line 167 | 3.086e18 |
    | Minimum wavelength for the simulated SED in meters | Line 168 | 9.2e-8 |
    | Maximum wavelength for the simulated SED in meters | Line 169 | 7.0e-7 |
    | Number of wavelength points in the simulated SED | Line 170 | 40 |

### Step 4: Run POLARIS using the Python scripts
Export two environmental variables:
* POLARIS_DIR: **absolute path** to the POLARIS repo, e.g. /home/POLARIS/
* TORCH_DIR: **absolute path** to the Torch repo, e.g. /home/Torch/

Then, run run\_polaris.py, using the arguments summarized in the table below. 

| Short option | Argument | Required | Default | Description |
|---|---|---|---|---|
| `-p` | `path_snapshot` | False | `'./'` | Absolute Path to the Torch snapshot file. |
| `-f` | `file_snapshot` | True if -s is not specified | `''` | Absolute Path to the file containing the list of snapshots to run POLARIS on. Each line should contain a comma-separated list of snapshot numbers (with no comma at the end). |
| `-l` | `line_snapshot` | False | `0` | The line number in the snapshot list file for POLARIS to process. |
| `-s` | `snapshot_number` | True if -f is not specified | `None` | The snapshot number to run POLARIS on. |
| `-o` | `output_directory`| True | `''` | The directory where POLARIS output will be saved. |
| `-ep` | `epsilon` | False | `1e-6` | The dust-to-gas ratio to use in POLARIS. |
| `-ns` | `nside` | False | `7` | The `nside` for the spherical array of detectors to use in POLARIS. |
| `-el` | `extra_label` | False | `'run1'` | The extra label to use in the POLARIS command file name. |
| `-sl` | `sed_label` | False | `None` | The SED label to use in the POLARIS command file name. |

**Important Notes**:
1. You can run either a single snapshot or a list of snapshots. One of the following two arguments (-s and -f) are required.
    * A single snapshot: using -s to specify the snapshot number (e.g. line 26 of /example_run/example_slurm_job.sh)
    * A list of snapshots: using -f to specify a comma-separated text file that contains the list of snapshots (e.g. line 28 of /example_run/example_slurm_job.sh)
2. /example_run/ contains an example shell script that shows how the arguments work, as well as an example file_snapshot. 
3. A quick test run without slurm would be to run lines 23 to 26 of /example_run/example_slurm_job.sh in the terminal. 
4. Most paths, including the environmental variables, have to be absolute.

