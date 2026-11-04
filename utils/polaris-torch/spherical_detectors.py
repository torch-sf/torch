import numpy as np
import healpy as hp
import os


def theta_phi_to_xy_angles(theta, phi):
    # Input: theta, phi in radians
    # Output: alpha, beta (rotation angles for POLARIS detectors) in radians
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    alpha_x = -np.arcsin(y)          # rotation about x [default 1st rotation axis by POLARIS]
    beta_y = np.arctan2(x,z)        # rotation about y [default 2nd rotation axis by POLARIS]

    alpha_x = alpha_x
    beta_y = beta_y

    return alpha_x, beta_y


def save_angle_detectors(dir_project, nside, save_coord=True, save_polaris_angle=True):
    # Save the coordinates / POLARIS detector angles for an array of detectors on spherical surface
    #
    # Input: nside is an integer
    #        dir_project is the directory in which the coordinates / POLARIS detectors angles are saved
    #        set save_coord to True to save spherical coordinates (theta, phi)
    #        set save_polaris_angle to True to save POLARIS detectors angles (alpha, beta)
    
    
    npix = hp.nside2npix(nside)  # Total number of pixels

    theta, phi = hp.pix2ang(nside, np.arange(npix))     # Get the (theta, phi) angles for each pixel center
    if save_coord:
        os.makedirs(dir_project+'/detector_angles/', exist_ok=True)
        np.savetxt(dir_project+f'/detector_angles/spherical_coordinates_detectors_nside{int(nside)}.txt',[theta, phi])

    alpha, beta = theta_phi_to_xy_angles(theta, phi)    # Convert spherical coordinates to 
    if save_polaris_angle:
        os.makedirs(dir_project+'/detector_angles/', exist_ok=True)
        np.savetxt(dir_project+f'/detector_angles/rotation_angles_detectors_nside{int(nside)}.txt',[alpha, beta])


    if not (save_coord or save_polaris_angle):
        return theta, phi, alpha, beta



def save_angle_detectors_lists(dir_project, nside_list, save_coord=True, save_polaris_angle=True):
    # Save the coordinates / POLARIS detector angles for a list of arrays of detectors on spherical surface
    # 
    # Input: dir_project is the directory in which the coordinates / POLARIS detectors angles are saved
    #        nside_list is a list
    #        set save_coord to True to save spherical coordinates (theta, phi)
    #        set save_polaris_angle to True to save POLARIS detectors angles (alpha, beta)

    os.makedirs(dir_project+'/detector_angles/', exist_ok=True) # Create the directory if it doesn't exist

    for nside in nside_list:
        npix = hp.nside2npix(nside)  # Total number of pixels
    
        theta, phi = hp.pix2ang(nside, np.arange(npix))     # Get the (theta, phi) angles for each pixel center
        if save_coord:
            np.savetxt(dir_project+f'/detector_angles/spherical_coordinates_detectors_nside{int(nside)}.txt',[theta, phi])
    
        alpha, beta = theta_phi_to_xy_angles(theta, phi)    # Convert spherical coordinates to 
        if save_polaris_angle:
            np.savetxt(dir_project+f'/detector_angles/rotation_angles_detectors_nside{int(nside)}.txt',[alpha, beta])
    
    
