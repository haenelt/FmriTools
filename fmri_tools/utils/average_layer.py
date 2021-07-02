# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys

# external inputs
import numpy as np
import nibabel as nb


def average_layer(img_input, path_output, basename_output, mode="mean"):
    """Average layer.
    
    This averages data across different layers or sessions. Input arrays should 
    be in mgh format. The output gets the suffix of the chosen mode.    

    Parameters
    ----------
    img_input : list
        List of input layers.
    path_output : str
        Path where output is saved.
    basename_output : str
        Basename of written output file.
    mode : str, optional
        Average mode (mean or median). The default is "mean".

    Returns
    -------
    None.
    
    """
    
    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)
    
    # initialise array
    data = nb.load(img_input[0])
    data_size = data.header["dims"][0]
    data_res = np.zeros((data_size, len(img_input)))
    
    # collect input arrays
    for i in range(len(img_input)):
        data_res[:, i] = nb.load(img_input[i]).get_fdata()[:, 0, 0]
    
    # average
    if mode == "mean":
        data_res = np.mean(data_res, axis=1)
    elif mode == "median":
        data_res = np.median(data_res, axis=1)
    else:
        sys.exit("Choose a valid mode!")
    
    # expand dimensions to match with input array
    data_res = np.expand_dims(data_res, axis=1)
    data_res = np.expand_dims(data_res, axis=1)
    
    # write output file
    output = nb.Nifti1Image(data_res, data.affine, data.header)
    nb.save(output, os.path.join(path_output, basename_output+"_"+mode+".mgh"))
