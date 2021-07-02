# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb


def get_laminar_profile(file_in, path_output, hemi, name_output, mode):
    """Get laminar profile.
    
    This function computes the laminar profile (mean, median, max, min) of data 
    sampled onto different cortical layers.    

    Parameters
    ----------
    file_in : list
        List of sampled data onto different cortical layers (overlay).
    path_output : str
        Path where to save laminar profile.
    hemi : str
        Hemisphere.
    name_output : str
        Output file name without file extension.
    mode : str
        mean, median, max, min.

    Returns
    -------
    None.
    
    """

    # make subfolders
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # get input header
    data_img = nb.load(file_in[0])

    # read input
    data = []
    for i in range(len(file_in)):
        data.append(np.squeeze(nb.load(file_in[i]).get_fdata()))

    # convert input into array
    data = np.array(data)

    # compute laminar profile
    if mode == "mean":
        data = np.mean(data, 0)
    elif mode == "median":
        data = np.median(data, 0)
    elif mode == "max":
        data = np.max(data, 0)
    elif mode == "min":
        data = np.min(data, 0)
    else:
        print("Choose a valid mode!")

    # expand axes to match input array
    data = np.expand_dims(data, 1)
    data = np.expand_dims(data, 1)

    # write output image
    output = nb.Nifti1Image(data, data_img.affine, data_img.header)
    nb.save(output,
            os.path.join(path_output, hemi + "." + name_output + ".mgh"))
