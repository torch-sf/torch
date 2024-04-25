#!/bin/sh

# Script that looks for last written checkpoint file in the output directory, finds the last part and plt files printed
# before this, and updates the flash.par file for restart. Keep in mind that this scripts overwrites the flash.par. 
# The script will not consider forced outputs, which are sometimes created when the code crashes. I suggest these cases 
# are handeled manually.
# Eric Andersson, 2024-March-21

# Prefix used for simulation 
filename_prefix=$(grep basenm flash.par | sed -n 's/.*"\([^"]*\)".*/\1/p')
output_dir=$(grep output_directory flash.par | sed -n 's/.*"\([^"]*\)".*/\1/p')
log_file=$(grep log_file flash.par | sed -n 's/.*"\([^"]*\)".*/\1/p')
echo "Filename prefix:" $filename_prefix
echo "Output directory:" $output_dir
echo "Logfile name:" $log_file

# Checkpoint file
chk_filename="$filename_prefix"hdf5_chk_0
last_chk=$(ls $output_dir | grep $chk_filename | tail -1 | sed "s|$chk_filename*\([0-9*]\)|\1|g")
echo "Last chk file number found: $last_chk"

# Find line in logfile where last file chk was created.
last_chk_linenumber=$(grep -n $chk_filename $log_file | cut -d: -f1 | tail -1)
last_chk_linenumber=$(($last_chk_linenumber + 6)) # Add 6 to capture plt and part files written at the same time as chk.

# Plot files
plt_filename="$filename_prefix"hdf5_plt_cnt_0
last_plt=$(head -n $last_chk_linenumber $log_file | grep $plt_filename | tail -1 | awk -F_ '{print $NF}' | sed 's/^0*//') 
echo "Last plt file number before last checkpoint: $last_plt"
next_plt=$(($last_plt +1 ))

# Particle files
part_filename="$filename_prefix"hdf5_part_0
last_part=$(head -n $last_chk_linenumber $log_file | grep $part_filename | tail -1 | awk -F_ '{print $NF}' | sed 's/^0*//') 
echo "Last part file number before last checkpoint: $last_part"
next_part=$(($last_part + 1))

# Note that this does not work on Mac (I did not check Windows). If used on Mac, add empty string after -i argument).
sed -i 's/restart *= *.false\./restart = .true./' flash.par
sed -i "s/\(checkpointFileNumber[[:space:]]*=[[:space:]]*\)[0-9]*/\1$last_chk/" flash.par
sed -i "s/\(plotFileNumber[[:space:]]*=[[:space:]]*\)[0-9]*/\1$next_plt/" flash.par
sed -i "s/\(particleFileNumber[[:space:]]*=[[:space:]]*\)[0-9]*/\1$next_part/" flash.par
