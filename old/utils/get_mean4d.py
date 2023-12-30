# -*- coding: utf-8 -*-

import os

import nibabel as nb
import numpy as np


def get_mean4d(file_in, path_output="", name_output="", write_output=False):
    """Get mean 4D.

    This function computes the mean time series of one or more time series.

    Parameters
    ----------
    file_in : list
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

    """

    # make subfolders
    if len(path_output) > 0:
        if not os.path.exists(path_output):
            os.makedirs(path_output)

    # get dimensions
    data_img = nb.load(file_in[0])

    res_array = np.zeros_like(data_img.get_fdata())
    for i in range(len(file_in)):
        res_array += nb.load(file_in[i]).get_fdata()

    res_array = res_array / len(file_in)

    # write mean time series
    output = nb.Nifti1Image(res_array, data_img.affine, data_img.header)
    if write_output:
        nb.save(output, os.path.join(path_output, "mean_" + name_output + ".nii"))

    return output
