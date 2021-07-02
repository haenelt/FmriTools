# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb


def get_series(file_in, path_out, name_output):
    """Get series.
    
    This function creates a 4D nifti time series from a set of 3D nifti files.

    Parameters
    ----------
    file_in : list
        Array of filenames containing single 3D nifti volumes.
    path_out : str
        Path where output is saved.
    name_output : str
        Basename of output 4D nifti file.

    Returns
    -------
    None.
    
    """

    # load 3D nifti to get array size
    data = nb.load(os.path.join(file_in[0]))
    data.header["dim"][0] = 4
    data.header["dim"][4] = len(file_in)

    res = np.zeros(data.header["dim"][1:5])
    for i in range(len(file_in)):
        img = nb.load(os.path.join(file_in[i])).get_fdata()
        res[:, :, :, i] = img

    output = nb.Nifti1Image(res, data.affine, data.header)
    nb.save(output, os.path.join(path_out, name_output + ".nii"))
