# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb
from nibabel.freesurfer.io import read_morph_data
from scipy.ndimage.filters import gaussian_filter


def map2grid(file_grid, file_input, sigma, path_output="", basename_output="", 
             binary=False, overwrite=True):
    """ Map to grid

    This script allows you to sample indexed morphological data onto the regular 
    grid. Optional, a gaussian filter can be applied to the output image.    

    Parameters
    ----------
    file_grid : str
        Filename of grid coordinate mapping.
    file_input : str
        Filename of morphological data or *.mgh data.
    sigma : float
        Standard deviation of Gaussian kernel.
    path_output : str, optional
        Path where output is saved. The default is "".
    basename_output : str, optional
        Basename of written output file. The default is "".
    binary : bool, optional
        Threshold output grid (for curvature file). The default is False.
    overwrite : bool, optional
        Write output to file. The default is True.

    Returns
    -------
    grid_array : ndarray
        File mapped onto array.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 01-11-2018             
    Last modified: 12-10-2020

    """
    
    # load data
    grid_img = nb.load(file_grid)
    grid_array = grid_img.get_fdata()
    if os.path.splitext(file_input)[1] == ".mgh":
        morph = nb.load(file_input).get_fdata()
    else:
        morph = read_morph_data(file_input)

    # sample data onto grid
    for i in range(np.size(grid_array,0)):
        for j in range(np.size(grid_array,1)):
            if grid_array[i,j] != 0:
                grid_array[i,j] = morph[grid_array[i,j].astype(int)]
    
    # gaussian filter (opt)
    if sigma != 0:
        order = 0
        mode = "reflect"
        truncate = 4.0
        grid_array = gaussian_filter(grid_array,
                                     sigma=sigma,
                                     order=order,
                                     mode=mode,
                                     truncate=truncate)
    
    # binary mode (opt)
    if binary is True:
        grid_array[grid_array > 0] = 1
        grid_array[grid_array != 1] = -1
    
    # write output data
    if overwrite:
        
        # make output folder
        if not os.path.exists(path_output):
            os.mkdir(path_output)
        
        if sigma == 0 and binary is True:
            filenameOUT = os.path.join(path_output,basename_output+"_grid_binary.nii")
        elif sigma == 0 and binary is False:
            filenameOUT = os.path.join(path_output,basename_output+"_grid.nii")
        elif sigma != 0 and binary is True:
            filenameOUT = os.path.join(path_output,basename_output+"_sigma"+str(sigma)+"_grid_binary.nii")
        else:
            filenameOUT = os.path.join(path_output,basename_output+"_sigma"+str(sigma)+"_grid.nii")
        
        output = nb.Nifti1Image(grid_array, grid_img.affine, grid_img.header)
        nb.save(output,filenameOUT)

    return grid_array
