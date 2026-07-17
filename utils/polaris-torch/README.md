## Utility to run Polaris on a Torch simulation snapshot
First written by Stefan Reissl (flash-to-polaris.py) and Brooke Polak (Polaris-torch.py)

Edited by Donglin Wu (all files)

To use:

1. Download and install Polaris: https://github.com/polaris-MCRT/POLARIS
2. Install python package healpy (and other more common dependencies) 
3. Edit script generate\_polaris\_scripts.py. This python file includes a function to produce command files for POLARIS. See comments in the script for more details.
4. Edit script run\_polaris.py. This python file uses generate\_polaris\_scripts.py to write command files for POLARIS and runs terminal command of POLARIS using subprocess. See comments in the script for more details.
Important parts to be edited: the path to the POLARIS directory, and the path to this directory.
5. Run run\_polaris.py, with the following arguments:

| Short option | Argument | Default | Description |
|---|---|---|---|
| `-p` | `path_snapshot` | `'./'` | Path to the Torch snapshot file. |
| `-f` | `file_snapshot` | `''` | Path to the file containing the list of snapshots to run POLARIS on. Each line should contain a comma-separated list of snapshot numbers. |
| `-l` | `line_snapshot` | `0` | The line number in the snapshot list file for POLARIS to process. |
| `-s` | `snapshot_number` | `None` | The snapshot number to run POLARIS on. |
| `-o` | `output_directory` | `''` | The directory where POLARIS output will be saved. |
| `-ep` | `epsilon` | `1e-6` | The dust-to-gas ratio to use in POLARIS. |
| `-ns` | `nside` | `7` | The `nside` for the spherical array of detectors to use in POLARIS. |
| `-el` | `extra_label` | `'run1'` | The extra label to use in the POLARIS command file name. |
| `-sl` | `sed_label` | `None` | The SED label to use in the POLARIS command file name. |
6. Have fun with the data.