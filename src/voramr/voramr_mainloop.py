from amuse.lab import *

from torch_amuse_flash.interface import Flash
from voramr.kdtree import (
    read_hdf5,
    build_kdtree,
    pickle_tree,
    unpickle_tree,
    interp_data,
    )
from voramr.hdf5_convert import (
    extract_data,
    rescale_coords_vels,
    write_corrected_file,
    )
from voramr.voramr_stdout import vprint

import numpy as np
from time import time

def get_ntasks_from_run_script(name="run.sh"):
    """formally -n is --ntasks, de facto same as nprocs"""
    # Unused in favor of same func in torch_mainloop.py - SCL
    n = None
    with open(name) as f:
        for line in f:
            w = line.split()
            if len(w) >= 3 and w[0] == '#SBATCH' and w[1] == '-n':
                assert n is None  # throw error if #SBATCH -n occurs >1x
                n = int(w[2])
    assert n is not None
    return n

def initialize_workers():
    # Unused in favor of same func in torch_mainloop.py - SCL
    vprint("Got ntasks from sbatch file: {}".format(USER['num_hy_workers']+1))
    vprint("Number of FLASH workers: {}".format(USER['num_hy_workers']))
    vprint("Initializing Hydro code...")
    # Converter for the N-body code.
    convert = nbody.nbody_to_si(1.0|units.parsec, 1000.0|units.MSun)
    # Converter for the hydro code.
    convert2 = generic_unit_converter.ConvertBetweenGenericAndSiUnits(1.0|units.cm, 1.0|units.g, 1|units.s)

    hydro = Flash(
        unit_converter=convert2,
        number_of_workers=USER['num_hy_workers'],
        redirection='file',
        redirect_stdout_file='voramr_worker.out',
        redirect_stderr_file='voramr_worker.err',
        )

    hydro.initialize_code()
    vprint("Hydro code initialized.")
    return hydro

def get_leaf_blocks(hydro, cellsPerBlock=16, numBlocks=None):
    """
    Acquires all FLASH blocks representing the computational domain.
    
    Arguments:
    hydro         - instance of AMUSE flash_worker.

    cellsPerBlock - number of cells in each direction of a FLASH block.

    numBlocks     - number of blocks in the computational domain. 
                This is passed in by the user, as long as numBlocks 
                is > the actual number of blocks in the domain, 
                this routine runs fine. 

    Returns:
    leaf_grids[:numblks] - list of active blocks. The slicing removes
                           any extra buffers that are present from
                           passing too many numBlocks.
    """
    vprint("Getting block data...")
    lim=cellsPerBlock
    lim3 = lim**3

    vprint("getting leaf_indices")
    all_grids=numBlocks # hard coded from torch_user.py. Can be set arbitratily large to
                        # accomodate any sized grid.
    [leaf_grids, block_array, num_leafs]= hydro.get_leaf_indices(list(range(all_grids)))
    numblks=num_leafs[0]
    
    leaf_grids = np.resize(leaf_grids,numblks*lim3)
    
    return leaf_grids[:numblks], block_array

def interpolate_fields(hydro, leaf_grids, proc_blocks, kdtree, nprocs, cellsPerBlock=16):
    """
    This subroutine loops over each block in the computational domain
    and does three things. 
    1.) The x,y,z coordinates of all cells
    within the block are extracted and transformed into a 3D numpy
    meshgrid. 
    2.) The coordinate mesh is passed to the kdtree 
    interpolator and nearest neighbor interpolation is performed on
    each cell within the mesh resulting in a 4D NxNxNxM matrix where
    N is the dimensionality of the FLASH block and M is the number of
    field values interpolated. 
    3.) The inerpolated field values
    are unraveled one by one from the 4D matrix and flattened into an
    ordered 1D array and fed back into FLASH. FLASH is smart enough
    to fill the 3D block matrix from a 1D array.
    
    Arguments:
    hydro      - instance of AMUSE flash_worker.
    
    leaf_grids - list of active blockIDs.

    proc_blocks = array of length number of processors, with value being number of blocks on that proc

    kdtree     - 3D tree object built previously from the input data, 
                 allows for nearest neighbor interpolation.
    """
    lim = cellsPerBlock
    a = np.empty(lim)
    a.fill(1)
    b = np.empty_like(a)
    b.fill(2)
    c = np.empty_like(a)
    c.fill(3)
    vprint("Getting {} block cell coords, interpolating, pass back to FLASH.".format(leaf_grids[-1]))
    # loop over all processors
    disp = 0
    procID = 0
    for i in range(nprocs): 
        # number of leaf blocks on this processor
        num_leaf_blks = proc_blocks[i]
        # loop over each blockID on each processor
        for j in range(num_leaf_blks): # Cycle over BlkIDs
            curr_blk_id = leaf_grids[j+disp]
            # Get x, y, z coordinates of cells in BlkID==leaf and myProc==leaf_proc
            x = np.array(hydro.get_1blk_cell_coords(a,curr_blk_id,procID,lim).value_in(units.cm))
            y = np.array(hydro.get_1blk_cell_coords(b,curr_blk_id,procID,lim).value_in(units.cm))
            z = np.array(hydro.get_1blk_cell_coords(c,curr_blk_id,procID,lim).value_in(units.cm))

            # Mesh coordinates into single 3D matrix object
            coords_mesh = np.meshgrid(x,y,z,indexing='ij')
            
            # Pass coordinate mesh into kdtree interpolation, get field values for each coord point.
            interp = interp_data(kdtree, coords_mesh)
            
            # Using interpolated data. Flatten to 1D to feed to Fortran.
            rho = interp[:,:,:,0].flatten(order='F') | units.g/units.cm**3
            eint = interp[:,:,:,1].flatten(order='F') | (units.cm**2)/(units.s**2)
            vx = interp[:,:,:,2].flatten(order='F') | units.cm/units.s
            vy = interp[:,:,:,3].flatten(order='F')	| units.cm/units.s
            vz = interp[:,:,:,4].flatten(order='F') | units.cm/units.s
            gpot = interp[:,:,:,5].flatten(order='F') | (units.cm**2)/(units.s**2)
            
            dataSize = cellsPerBlock
            
            # Feed field data to FLASH.
            #Fortran can properly populate its NxNxN matrices when given a 1xN^3 matrix
            hydro.set_block_state(curr_blk_id, procID, dataSize, rho, vx, vy, vz, eint, gpot)
        disp += num_leaf_blks
        procID += 1
    #vprint("Done setting blocks. Total blocks: ", leaf)
    

def run_flash(user_initial_conditions, user_parameters):
    """
    """
    # This is not used now that VorAMR is embbedded within Torch. Instead,
    # the run_flash() function defined in torch_mainloop.py is used. - SCL
    global USER
    USER = user_parameters()

    if(USER['convert_file']):
        coords, vels, dens, mass, eint, gpot = extract_data(USER['source_file'],
                                                      apply_consts=True)
        coords_cor, vels_cor = rescale_coords_vels(coords, vels, mass, use_com_coords=False)
        write_corrected_file(USER['input_file'], coords_cor, vels_cor, dens, mass, eint, gpot)
        
        coords, field_set = read_hdf5(USER['input_file'])
    else:
        coords, field_set = read_hdf5(USER['source_file'])

    vprint("Building field interpolator.")
    kdtree = build_kdtree(coords, field_set)
    if(USER['pickle_kdtree']):
        pickle_tree(kdtree, USER['pickle_file_name'])
        
    vprint("Running Flash.")
    hydro = initialize_workers()
    leaf_blocks = get_leaf_blocks(hydro, cellsPerBlock=USER['cellsPerBlock'], numBlocks=USER['numBlocks'])

    interpolate_fields(hydro, leaf_blocks, kdtree, cellsPerBlock=USER['cellsPerBlock'])

    vprint("Trying write_chpt()")
    hydro.write_chpt()
    vprint("Trying IO_out('chk', 5)")
    hydro.IO_out('chk', 5)
    try:
        vprint("Evolving Flash...")
        evolve(hydro)
        
    finally:
        pass
