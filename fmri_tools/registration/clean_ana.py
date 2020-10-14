# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb


def clean_ana(input, min_value, new_range, overwrite=True):
    """ Clean ana

    This function removes ceiling values from a computed T1 map of an mp2rage 
    acquisition. Low intensity values are removed and the data range is 
    normalised to a defined new range. The input file should be either a nifti 
    or a compressed nifti file.    

    Parameters
    ----------
    input : str
        Filename of input image.
    min_value : float
        Threshold of low intensity values.
    new_range : float
        Arbitrary new data range.
    overwrite : bool, optional
        If set, the input image is overwritten with the cleaned data set. The 
        default is True.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 20-10-2019        
    Last modified: 12-10-2020
    
    """
    
    # load data
    data = nb.load(input)
    data_array = data.get_fdata()
    
    # remove ceiling
    data_array[data_array == 0] = np.max(data_array)
    
    # remove low intensity values
    data_array[data_array <= min_value] = 0
    data_array = data_array - np.min(data_array[data_array != 0])
    data_array[data_array <= 0] = 0
    
    # normalise to new data range
    data_array = data_array / np.max(data_array) * new_range
    
    # output cleaned dataset
    output = nb.Nifti1Image(data_array, data.affine, data.header)
    if overwrite:
        nb.save(output, input)
    else:
        path = os.path.dirname(input)
        if os.path.splitext(input)[1] == ".gz":
            basename = os.path.splitext(os.path.splitext(os.path.basename(input))[0])[0]
            nb.save(output, os.path.join(path, basename+"_clean.nii.gz"))
        else:
            basename = os.path.splitext(os.path.basename(input))[0]
            nb.save(output, os.path.join(path, basename+"_clean.nii"))
            