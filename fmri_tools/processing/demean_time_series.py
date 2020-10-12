# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb


def demean_time_series(img_input, path_output="", name_output="", 
                       write_output=False):
    """
    This function demeans each voxel time series. Input is either a 4d nifti or 
    compressed nifti file.
    Inputs:
        *img_input: 4d nifti volume or string to filename.
        *path_output: path where output is saved.
        *name_output: basename of output.
        *write_output: write nifti volume.
    Outputs:
        *output: demeaned 4d nifti volume.
    
    created by Daniel Haenelt
    Date created: 24-10-2019
    Last modified: 12-10-2020
    """

    # load data
    if isinstance(img_input, nb.Nifti1Image):
        data_array = img_input.get_fdata()
    elif isinstance(img_input, str):
        img_input = nb.load(img_input)
        data_array = img_input.get_fdata()
    else:
        print("Input must be either string or instance of nibabel class")
        return
    
    # get mean of each voxel time series
    data_mean = np.mean(data_array, axis=3)
    
    # demean time series
    for i in range(np.shape(data_array)[3]):
        data_array[:,:,:,i] = ( data_array[:,:,:,i] - data_mean )/ data_mean * 100
       
    # write output    
    output = nb.Nifti1Image(data_array, img_input.affine, img_input.header)
    if write_output:
        nb.save(output, os.path.join(path_output, "demean_"+name_output+".nii"))
    
    return output
