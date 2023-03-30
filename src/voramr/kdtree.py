from scipy.interpolate import NearestNDInterpolator
import h5py
import numpy as np
import pickle

def read_hdf5(file_path):
    """
    Opens an hdf5 data file and extracts the coordinates
    of the Voronoi cells (or other data structures)
    as well as the associated gas desnity, 
    internal energy, and velocity field values.
    
    This function is written specifically to expect AREPO output
    data. Other data structures will need to have tailored calls
    to extract the data from the hdf5 file structure. Hopefully
    you'll only have to edit this file though!

    Arguments:
    file_path - input AREPO hdf5 file path

    Returns:
    coords    - list containing 3 1xN numpy arrays representing the 
                x, y, z coordinate sets for all Voronoi cells.
    field_set - 5xN numpy array where the columns are the separate
                field values of interest (density, internal energy,
                velx, vely, velz).
    """
    f = h5py.File(file_path,'r')
    coords_set = np.array(f["PartType0"]["Coordinates"])
    x = coords_set[:,0]
    y = coords_set[:,1]
    z = coords_set[:,2]
    coords = [x,y,z]

    # Extract field values from AREPO HDF5
    density_set = np.array(f["PartType0"]["Density"])
    intEner_set = np.array(f["PartType0"]["InternalEnergy"])
    mass_set = np.array(f["PartType0"]["Masses"])
    velocity_set = np.array(f["PartType0"]["Velocities"])
    # We need to separate velocity vector, need to interpolate each component
    velx = velocity_set[:,0]
    vely = velocity_set[:,1]
    velz = velocity_set[:,2]
    gpot = np.array(f["PartType0"]["Potential"]) 

    field_set = np.stack((density_set, intEner_set, velx, vely, velz, gpot), axis=-1)

    return coords, field_set

def build_kdtree(coords, field_set):
    """
    Builds KDtree object from N coordinates and m x N field values.
    Each leaf of the tree corresponds to a Voronoi cell center with m 
    field values associated with that cell.

    Arguments:
    coords    - coordinate set [x, y, z] where x,y,z are 1 x N numpy arrays.
    field_set - m x N array of field values.

    Returns:
    tree      - kdtree object with field value matrix on each leaf node.
    """
    # Create interpolator object (tree)
    tree = NearestNDInterpolator(list(zip(coords[0], coords[1], coords[2])), field_set)
    return tree


def pickle_tree(tree, file_out_name):
    """
    Pickles any object. Syntax tailored specifically for tree
    structures produced by the functions in voramr_kdtree.py.
    The pickled tree object can be accessed an interpolated from
    without having to re-generate the tree.

    Arguments:
    tree          - tree (or any other) object to be pickled.
    file_out_name - file name to be saved in current directory.
    """
    # Pickling tree
    file_w = open(file_out_name, 'wb')
    pickle.dump(tree, file_w)
    file_w.close()
    return

def unpickle_tree(file_in_name):
    """
    Unpickles an object and returns it.

    Arguments:
    file_in_name - name/path of pickle file.

    Returns:    
    tree_struct - unpickled tree (or any other object).
    """
    file_r = open(file_in_name, 'rb')
    tree_struct = pickle.load(file_r)
    file_r.close()
    return tree_struct


def interp_data(tree_struct, cell_coords):
    """
    Performs Nearest Neighbor interpolation between single cell
    3d coordinates and a KDtree and returns the values associated
    with the nearest tree leaf (Voronoi cell center).

    Arguments:
    tree_struct - KDtree structure.
    cell_coords - [x, y, z] coordinates of single cell.

    Returns:
    interp_data - field values of nearest neighbor between cell 
                  coordinates and a KDtree leaf.
    """
    interp_data = tree_struct(cell_coords[0], cell_coords[1], cell_coords[2])
    return interp_data
