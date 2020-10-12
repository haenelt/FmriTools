# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb

# local inputs
from fmri_tools.io.get_filename import get_filename


def crop_coordinate_mapping(input, pad=0, overwrite_file=True, path_output=""):
    """
    Crops a padded coordinate mapping. The output file can either overwrite the 
    input file or a new file is created with a suffix in a defined output 
    directory.
    Inputs:
        *input: input file.
        *pad: image padding size.
        *overwrite_file: output file overwrites input file.
        *path_output: path where output is saved if input file is not overwritten.

    created by Daniel Haenelt
    Date created: 21-11-2018             
    Last modified: 12-10-2020
    """

    # define output folder
    if path_output is not None:
        if not os.path.exists(path_output):
            os.makedirs(path_output)

    # get input path and file name
    path, file, ext = get_filename(input)

    # load data
    data_img = nb.load(input)
    data_array = data_img.get_fdata()

    # get matrix size
    x_size = np.size(data_array,0)
    y_size = np.size(data_array,1)
    z_size = np.size(data_array,2)

    # crop image matrix
    data_array = data_array[pad:x_size-pad,pad:y_size-pad,pad:z_size-pad,:]

    # write cropped coordinate mapping
    output = nb.Nifti1Image(data_array, data_img.affine, data_img.header)
    output.set_data_dtype(np.float)

    # write coordinate mapping for each time point
    if overwrite_file is True:
        os.remove(input)
        nb.save(output, input)
    else:
        fileOUT = os.path.join(path_output, file+'_crop'+ext)
        nb.save(output,fileOUT)
