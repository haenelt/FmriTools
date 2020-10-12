# -*- coding: utf-8 -*-

# python standard library inputs
import os
import copy

# external inputs
import numpy as np
import nibabel as nb
    
    
def get_std(input, path_output, name_output, set_outlier=None):
    """
    This function computes the standard deviation of one or more time series.
    Inputs:
        *input: single file or list of files.
        *path_output: path where to save mean image
        *name_output: output file name without file extension.
        *set_outlier: can be nan, zero or None.
        
    created by Daniel Haenelt
    Date created: 04-02-2019         
    Last modified: 12-10-2020
    """
    
    # make subfolders
    if not os.path.exists(path_output):
        os.makedirs(path_output)
    
    if len(np.shape(input)) > 0:
        
        # get dimensions
        data_img = nb.load(input[0])
        
        x_size = data_img.header["dim"][1]
        y_size = data_img.header["dim"][2]
        z_size = data_img.header["dim"][3]
        
        t_size = 0
        for i in range(len(input)):
            data_img = nb.load(input[i])
            t_size += data_img.header["dim"][4]
        
        data_all_array = np.zeros([x_size, y_size, z_size, t_size])
        t_start = 0
        for i in range(len(input)):
            data_img = nb.load(input[i])
            data_array = data_img.get_fdata()
            t_end = data_img.header["dim"][4] + t_start                  
            data_all_array[:,:,:,t_start:t_end] = data_array
            t_start = copy.deepcopy(t_end)            
                
    else:
               
        # load data
        data_img = nb.load(input)
        data_all_array = data_img.get_fdata()
    
    # calculate std
    data_std_array = np.std(data_all_array,axis=3)
    
    if set_outlier == "nan":
        data_std_array[data_std_array == 0] = np.nan # set zeroes to nan
    elif set_outlier == "zero":
        data_std_array[data_std_array == 0] = 0
       
    # write output
    if len(np.shape(input)) > 0:
        data_img = nb.load(input[0])
    else:
        data_img = nb.load(input)
    data_img.header["dim"][0] = 3
    data_img.header["dim"][4] = 1
    
    # write mean image
    mean_img = nb.Nifti1Image(data_std_array, data_img.affine, data_img.header)
    nb.save(mean_img, os.path.join(path_output,"std_"+name_output+".nii"))
