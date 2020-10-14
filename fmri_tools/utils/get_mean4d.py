# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb

   
def get_mean4d(input, path_output="", name_output="", write_output=False):
    """ Get mean 4D

    This function computes the mean time series of one or more time series.    

    Parameters
    ----------
    input : list
        List of 4d nifti files.
    path_output : str, optional
        Path where to save mean image The default is "".
    name_output : str, optional
        Output file name without file extension. The default is "".
    write_output : bool, optional
        Write nifti volume. The default is False.

    Returns
    -------
    output : niimg
        Mean time series.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 31-10-2019         
    Last modified: 12-10-2020
    
    """
    
    # make subfolders
    if len(path_output) > 0:
        if not os.path.exists(path_output):
            os.makedirs(path_output)
        
    # get dimensions
    data_img = nb.load(input[0])
    
    res_array = np.zeros_like(data_img.get_fdata())
    for i in range(len(input)):
        res_array += nb.load(input[i]).get_fdata()

    res_array = res_array / len(input)
    
    # write mean time series
    output = nb.Nifti1Image(res_array, data_img.affine, data_img.header)
    if write_output:
        nb.save(output,os.path.join(path_output,"mean_"+name_output+".nii"))
    
    return output
