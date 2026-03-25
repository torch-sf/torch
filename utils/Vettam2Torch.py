# Written by Shyam Menon
# This script converts vanilla vettam to Torch-compatible vettam 
# Variable namees are changed
# NOTE: this script must be updated, see diff files between torch-vettam and vettam
# attached to the vettam pull request for all changes. 

import os
import re
import tarfile
import shutil


# Define the replacements
replacements = {
    r"\bEUVRATE_PART_PROP\b": "NION_PART_PROP",
    r"\bTHERMPION_PART_PROP\b": "EION_PART_PROP",
    r"\bLUMFUV_PART_PROP\b": "NPEP_PART_PROP",
    r"\bIH_VAR\b": "PHIO_VAR",
    r"\bIONH_VAR\b": "PHHE_VAR",
    r"\bHP_SPEC\b": "IHP_SPEC",
    r"\bH_SPEC\b": "IHA_SPEC",
    r"\bUPED_VAR\b": "PEFL_VAR"
}

def process_file(filepath):
    """Read a file, replace specified strings while ensuring whole-word matches, and save changes."""
    with open(filepath, "r") as file:
        content = file.read()
    
    modified = False
    for old, new in replacements.items():
        if re.search(old, content):
            content = re.sub(old, new, content)
            modified = True
    
    if modified:
        with open(filepath, "w") as file:
            file.write(content)
        print(f"Updated: {filepath}")

def scan_and_replace(directory):
    """Recursively scan a directory and process all .f90 files in all subdirectories."""
    for root, _, files in os.walk(directory):
        for file in files:
            if file.lower().endswith(".f90"):
                process_file(os.path.join(root, file))

def PEFL_def_change():
    """
    Modify the definition of PEFL_VAR in rt_dustTerms.F90. Default VETTAM stores the radiation energy density absorbed.
    To convert this to a flux in Habing units, we need to multiply by cl_speedlt/(4*PI)/2.1e-4.
    """
    filename = "rt_dustTerms.F90"

    # Read the file
    with open(filename, "r") as f:
        lines = f.readlines()

    # Define patterns and replacements
    pattern_if = re.compile(r"if\(current_band\s*\.eq\.\s*'PE'\)")
    replacement_if = "if(current_band .eq. 'FUV')"

    pattern_soln = re.compile(
        r"(solnData\(PEFL_VAR,i,j,k\)\s*=\s*solnData\(PEFL_VAR,i,j,k\)\s*\+.*?dt)"
    )
    replacement_soln = r"\1 * rt_speedlt/1.6e-3"

    # Modify the lines
    modified_lines = []
    inside_block = False

    for line in lines:
        if "#ifdef PEFL_VAR" in line:
            inside_block = True

        if inside_block:
            line = pattern_if.sub(replacement_if, line)
            line = pattern_soln.sub(replacement_soln, line)

        if "#endif" in line:
            inside_block = False

        modified_lines.append(line)

    # Write the modified content back to the file
    with open(filename, "w") as f:
        f.writelines(modified_lines)

    print(f"Modifications applied to {filename}")

def NPEP_def_change():
    file_path = "rt_sinkInject.F90"
    # Open the file in read mode
    with open(file_path, 'r') as file:
        content = file.read()

    # Define the block to be modified
    block_pattern = r'(#elif NION_PART_PROP.*?)(!Star case.*?#else)'  # Regex to capture the block

    # Define the replacement for the block
    replacement = '''#elif NION_PART_PROP
    integer, parameter :: gather_nprops = 3
    integer, dimension(gather_nprops), save :: gather_propinds = &
    (/ integer :: NION_PART_PROP, NPEP_PART_PROP, EPEP_PART_PROP /)
  !Star case (e.g.: /scratch/ek9/sm5890/flash_newnew/flash-rsaa/source/Particles/ParticlesMain/Sink/StellarEvolution)
#else'''

    # Replace the block using regex
    modified_content = re.sub(block_pattern, replacement, content, flags=re.DOTALL)

    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.write(modified_content)    

    # Open the file in read mode
    with open(file_path, 'r') as file:
        content = file.read()

    # Replace all occurrences of the lum_sink line throughout the entire file
    modified_content = re.sub(r'lum_sink = particles_global\(NPEP_PART_PROP,p\)', 
                     'lum_sink = particles_global(NPEP_PART_PROP,p)*particles_global(EPEP_PART_PROP,p)', 
                     content)

    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.write(modified_content)    

    print(f"Modifications applied to {file_path}")

if __name__ == "__main__":
    #Change directory to the RadTrans/VETTAM directory
    os.chdir("flash/source/physics/RadTrans/RadTransMain/VETTAM")
    repo_directory = os.getcwd()  # User-specified directory
    #Make variable substitutions
    scan_and_replace(repo_directory)
    #Make PEFL_VAR definition changes
    PEFL_def_change()
    #Make sure definition of the FUV band luminosity is appropriate for torch (LUMFUV -> EPEP*NPEP)
    NPEP_def_change()

    print("Converting VETTAM to Torch complete.")
