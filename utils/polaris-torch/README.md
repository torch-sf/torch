## Utility to run Polaris on a Torch simulation snapshot
First written by Stefan Reissl (flash-to-polaris.py) and Brooke Polak (Polaris-torch.py, generate\_polaris\_scripts.py)

Edited by Donglin Wu

To use:

1. Download and install Polaris: https://github.com/polaris-MCRT/POLARIS
2. Install python package healpy (and other more common dependencies) 
3. Edit script generate\_polaris\_scripts.py. This python file includes a function to produce command files for POLARIS. See comments in the script for more details.
4. Edit script run\_polaris.py. This python file uses generate\_polaris\_scripts.py to write command files for POLARIS and runs terminal command of POLARIS using subprocess. See comments in the script for more details.
Important lines to be edited: the path to the POLARIS directory (dir_polaris), and the path to this directory (dir_run).
5. Run run\_polaris.py, with the arguments in the table below.
6. Have fun with the data.


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

**Note**:
<ol type="A">
  <li> You can run either a single snapshot (using -s) or a list of snapshots (using -f and specifying a text file that contains the list). One of the two arguments are required to be specified. </li>
  <li> The paths have to be absolute. </li>
  <li> /example_run/ contains an example shell script that shows how the arguments work, as well as an example file_snapshot.  </li>
</ol>

