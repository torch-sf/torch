## Utility to run Polaris on a Torch simulation snapshot
Written by Stefan Reissl (flash-to-polaris.py) and Brooke Polak (Polaris-torch.py)

To use:

1. Download and install Polaris: https://github.com/polaris-MCRT/POLARIS
2. Copy a torch snapshot to desired run directory, along with these three python files
3. Edit example script generate\_polaris\_sim.py.
4. Run python generate\_polaris\_sim.py. This will produce the file torch\_cmd.
5. To run the polaris sim, execute polaris torch\_cmd
6. Have fun with the data (hopefully)!
