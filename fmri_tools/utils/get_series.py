# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb


def get_series(input, path_out, name_output):
    """
    This function creates a 4D nifti time series from a set of 3D nifti files.
    Inputs:
        *input: array of filenames containing single 3D nifti volumes.
        *path_output: path where output is saved.
        *name_output: basename of output 4D nifti file.
        
    created by Daniel Haenelt
    Date created: 28-06-2020       
    Last modified: 12-10-2020
    """

    # load 3D nifti to get array size
    data = nb.load(os.path.join(input[0]))
    data.header["dim"][0] = 4
    data.header["dim"][4] = len(input)

    res = np.zeros(data.header["dim"][1:5])
    for i in range(len(input)):

        img = nb.load(os.path.join(input[i])).get_fdata() 
        res[:,:,:,i] = img
    
    output = nb.Nifti1Image(res, data.affine, data.header)
    nb.save(output,os.path.join(path_out,name_output+".nii"))
    