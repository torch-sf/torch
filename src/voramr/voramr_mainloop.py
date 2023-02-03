from amuse.lab import *
#from amuse.datamodel import Particles

from amuse.community.voramr.interface import Flash
from voramr_kdtree import (
    read_arepo_hdf5,
    build_kdtree,
    pickle_tree,
    unpickle_tree,
    interp_data,
    )
from voramr_convert import (
    extract_data,
    rescale_coords_vels,
    write_corrected_file,
    )
from voramr_stdout import vprint

import numpy as np
from time import time

def get_ntasks_from_run_script(name="run.sh"):
    """formally -n is --ntasks, de facto same as nprocs"""
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
    vprint("Getting block data...")
    lim=cellsPerBlock
    lim3 = lim**3
    # nprocs = USER['num_hy_workers']
    nprocs = hydro.get_number_of_procs()
    vprint("true nprocs:", nprocs)
    all_grids = hydro.get_number_of_grids(nprocs)
    vprint("true all grids:", all_grids)
    nprocs = 4 #get_number_of_procs()
    # get all blocks (grids)                                                                                  
    all_grids = hydro.get_number_of_grids(nprocs)
    # This returns 0 and I'm not sure why.
    vprint("all grids from nprocs=4:", all_grids)
    leaf_grids = np.zeros((all_grids*lim3))
    block_array = np.zeros((all_grids*lim3))

    vprint("getting leaf_indices")
    all_grids=numBlocks # hard code number of blocks testing SCL 08/31
    [leaf_grids, block_array, num_leafs]= hydro.get_leaf_indices(list(range(all_grids)))
    numblks=num_leafs[0]
    
    leaf_grids = np.resize(leaf_grids,numblks*lim3)
    
    return leaf_grids[:numblks]

def interpolate_fields(hydro, leaf_grids, kdtree, cellsPerBlock=16):
    lim = cellsPerBlock
    a = np.empty(lim)
    a.fill(1)
    b = np.empty_like(a)
    b.fill(2)
    c = np.empty_like(a)
    c.fill(3)
    vprint("Getting {} block cell coords, interpolating, pass back to FLASH.".format(leaf_grids[-1]))
    for leaf in leaf_grids: # Cycle over BlkIDs
        # Get x, y, z coordinates of cells in BlkID==leaf
        x = np.array(hydro.get_1blk_cell_coords(a,leaf,lim).value_in(units.cm))
        y = np.array(hydro.get_1blk_cell_coords(b,leaf,lim).value_in(units.cm))
        z = np.array(hydro.get_1blk_cell_coords(c,leaf,lim).value_in(units.cm))

        # Mesh coordinates into single 3D matrix object
        coords_mesh = np.meshgrid(x,y,z,indexing='ij')
        
        # Pass coordinate mesh into kdtree interpolation, get field values for each coord point.
        interp = interp_data(kdtree, coords_mesh)
        
        i = np.arange(1,17) # indicies 1-16 (Fortran indexes from 1)
        j, k = i.copy(), i.copy()
        index_mesh = np.meshgrid(i, j, k, indexing='xy')
        BlkIndex = leaf
        nproc = 1
        # --------
        
        # Using interpolated data. Flatten to 1D to feed to Fortran.
        rho = interp[:,:,:,0].flatten(order='F') | units.g/units.cm**3
        eint = interp[:,:,:,1].flatten(order='F') | (units.cm**2)/(units.s**2)
        vx = interp[:,:,:,2].flatten(order='F') | units.cm/units.s
        vy = interp[:,:,:,3].flatten(order='F')	| units.cm/units.s
        vz = interp[:,:,:,4].flatten(order='F') | units.cm/units.s
        gpot = interp[:,:,:,5].flatten(order='F') | (units.cm**2)/(units.s**2)
        
        dataSize = 16
        
        # Feed field data to FLASH.
        #Fortran can properly populate its NxNxN matrices when given a 1xN^3 matrix
        hydro.set_block_state(BlkIndex, dataSize, rho, vx, vy, vz, eint, gpot)
        
    vprint("Done setting blocks. Total blocks: ", leaf)
    

def run_flash(user_initial_conditions, user_parameters):
    """
    """
    global USER
    USER = user_parameters()

    if(USER['convert_file']):
        coords, vels, dens, mass, eint, gpot = extract_data(USER['source_file'],
                                                      apply_consts=True)
        coords_cor, vels_cor = rescale_coords_vels(coords, vels, mass, apply_consts=True, use_com_coords=False)
        write_corrected_file(USER['input_file'], coords_cor, vels_cor, dens, mass, eint, gpot)
        
        coords, field_set = read_arepo_hdf5(USER['input_file'])
    else:
        coords, field_set = read_arepo_hdf5(USER['source_file'])

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
