# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb

# local inputs
from fmri_tools.io import get_filename


def write_mgh(arr, affine, header, file_out):
    """ Write MGH

    This function adds two empty dimensions to an array and saves it as a
    freesurfer mgh surface file.        

    Parameters
    ----------
    arr : ndarray
        Image array.
    affine : ndarray
        Affine transformation matrix.
    header : MGHHeader
        Image header.
    file_out : str
        Filename of output file.

    Raises
    ------
    ValueError
        If `file_out` is not a string or has a file extension which is not 
        supported.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 25-08-2020
    Last modified: 20-10-2020

    """

    # check filename
    if isinstance(file_out, str):
        if not file_out.endswith("mgh"):
            raise ValueError("Currently supported file formats is mgh.")
    else:
        raise ValueError("Filename must be a string!")
    
    # make output folder
    path_output, _, _ = get_filename(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # add empty dimensions
    arr = np.expand_dims(arr, axis=1)
    arr = np.expand_dims(arr, axis=1)

    # write output
    output = nb.Nifti1Image(arr, affine, header)
    nb.save(output, file_out)
