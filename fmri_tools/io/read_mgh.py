# -*- coding: utf-8 -*-

# external inputs
import numpy as np
import nibabel as nb


def read_mgh(file_in):
    
    """ Read MGH
    
    This function reads an surface mgh file and removes empty dimensions to
    return an array.    

    Parameters
    ----------
    file_in : str
        Filename of input file.

    Returns
    -------
    arr : ndarray
        Image array.
    affine : ndarray
        Affine transformation matrix.
    header : MGHHeader
        Image header.
    
    Notes
    -------
    created by Daniel Haenelt
    Date created: 25-08-2020
    Last modified: 20-10-2020

    """
    
    # get header
    header = nb.load(file_in).header
    affine = nb.load(file_in).affine
    
    # get data
    arr = nb.load(file_in).get_fdata()
    arr = np.squeeze(arr)
    
    return arr, affine, header
