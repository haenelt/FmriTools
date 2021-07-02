# -*- coding: utf-8 -*-

# external inputs
import numpy as np
import nibabel as nb


def read_mgh(file_in):
    """Read MGH.

    This function reads a surface mgh file and removes empty dimensions from the
    data array.

    Parameters
    ----------
    file_in : str
        File name of input file.

    Raises
    ------
    ValueError
        If `file_in` is not a string or has a file extension which is not 
        supported.

    Returns
    -------
    arr : ndarray
        Image array.
    affine : ndarray
        Affine transformation matrix.
    header : MGHHeader
        Image header.

    """
    
    # check filename
    if isinstance(file_in, str):
        if not file_in.endswith("mgh"):
            raise ValueError("Currently supported file format is mgh.")
    else:
        raise ValueError("Filename must be a string!")
    
    # get header
    header = nb.load(file_in).header
    affine = nb.load(file_in).affine
    
    # get data
    arr = nb.load(file_in).get_fdata()
    arr = np.squeeze(arr)
    
    return arr, affine, header
