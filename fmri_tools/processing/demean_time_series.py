# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb


def demean_time_series(img_input, path_output="", name_output="",
                       write_output=False):
    """Demean time series.

    This function demeans each voxel time series. Input is either a 4d nifti or 
    compressed nifti file.    

    Parameters
    ----------
    img_input : niimg
        4d nifti volume or string to filename.
    path_output : str, optional
        Path where output is saved. The default is "".
    name_output : str, optional
        Basename of output. The default is "".
    write_output : bool, optional
        Write nifti volume. The default is False.

    Returns
    -------
    output : niimg
        Demeaned 4d nifti volume.
    
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
        data_array[:, :, :, i] = (data_array[:, :, :, i] - data_mean) / data_mean * 100

    # write output    
    output = nb.Nifti1Image(data_array, img_input.affine, img_input.header)
    if write_output:
        nb.save(output,
                os.path.join(path_output, "demean_" + name_output + ".nii"))

    return output
