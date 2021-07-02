# -*- coding: utf-8 -*-

# python standard library inputs
import os
import copy

# external inputs
import numpy as np
import nibabel as nb


def get_mean(file_in, path_output, name_output, method="mean"):
    """Get mean.

    This function computes the mean image of one or more time series.

    Parameters
    ----------
    file_in : str or list
        Single file or list of files.
    path_output : str
        Path where to save mean image
    name_output : str
        Output file name without file extension.
    method : str, optional
        Can be either mean or median. The default is "mean".

    Raises
    ------
    ValueError
        If `method` is invalid.

    Returns
    -------
    None.

    """

    # make subfolders
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    if len(np.shape(file_in)) > 0:

        # get dimensions
        data_img = nb.load(file_in[0])

        x_size = data_img.header["dim"][1]
        y_size = data_img.header["dim"][2]
        z_size = data_img.header["dim"][3]

        t_size = 0
        for i in range(len(file_in)):
            data_img = nb.load(file_in[i])
            t_size += data_img.header["dim"][4]

        data_all_array = np.zeros([x_size, y_size, z_size, t_size])
        t_start = 0
        for i in range(len(file_in)):
            data_img = nb.load(file_in[i])
            data_array = data_img.get_fdata()
            t_end = data_img.header["dim"][4] + t_start
            data_all_array[:, :, :, t_start:t_end] = data_array
            t_start = copy.deepcopy(t_end)

    else:

        # load data
        data_img = nb.load(file_in)
        data_all_array = data_img.get_fdata()

    # calculate mean
    if method == "mean":
        data_mean_array = np.mean(data_all_array, axis=3)
    elif method == "median":
        data_mean_array = np.median(data_all_array, axis=3)
    else:
        raise ValueError("Choose a valid mean type!")

    # write output
    if len(np.shape(file_in)) > 0:
        data_img = nb.load(file_in[0])
    else:
        data_img = nb.load(file_in)
    data_img.header["dim"][0] = 3
    data_img.header["dim"][4] = 1

    # write mean image
    mean_img = nb.Nifti1Image(data_mean_array, data_img.affine, data_img.header)
    nb.save(mean_img, os.path.join(path_output, "mean_" + name_output + ".nii"))
