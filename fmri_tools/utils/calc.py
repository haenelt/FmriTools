# -*- coding: utf-8 -*-
"""Different image manipulation tools."""

import os
import sys

import nibabel as nb
import numpy as np

from ..utils.bias import remove_bias_ants

__all__ = ["estimate_t1w", "mip", "remove_nans", "multiply_images", "average_layer"]


def estimate_t1w(file_bold, file_vaso, file_out, apply_bias):
    """Compute a T1-weighted image from nulled and not-nulled time series. Be
    aware that the computation takes all time points. Therefore, volumes not in
    steady-state should be removed before running this function.

    Parameters
    ----------
    file_bold : str
        File name of not-nulled time series.
    file_vaso : str
        File name of nulled time series.
    file_out : str
        File name of T1-weighted image.
    apply_bias : bool
        Apply a bias field correction.

    Returns
    -------
    None.

    """

    print("Compute T1-weighted image from")
    print("BOLD: " + file_bold)
    print("VASO: " + file_vaso)

    # load data
    bold = nb.load(file_bold)
    vaso = nb.load(file_vaso)
    arr_bold = bold.get_fdata()
    arr_vaso = vaso.get_fdata()

    bold_dim = bold.header["dim"]
    nx = bold_dim[1]
    ny = bold_dim[2]
    nz = bold_dim[3]
    nt = bold_dim[4]

    arr_combined = np.zeros((nx, ny, nz, 2 * nt))
    arr_combined[:, :, :, :nt] = arr_bold
    arr_combined[:, :, :, nt:] = arr_vaso

    arr_abs = np.abs(np.mean(arr_combined, axis=3))
    arr_std = np.std(arr_combined, axis=3)
    arr_std[arr_std == 0] = 1
    arr_std[~np.isfinite(arr_std)] = 1
    arr_res = arr_abs / arr_std

    # write output
    bold.header["dim"][0] = 3
    bold.header["dim"][4] = 1
    output = nb.Nifti1Image(arr_res, bold.affine, bold.header)
    nb.save(output, file_out)

    # bias field correction to epi
    if apply_bias:
        print("Apply bias field correction")
        remove_bias_ants(file_out, file_out, save_bias=False)


def mip(file_in, file_out, size, axis=2, mode="min"):
    """Minimum/Maximum intensity projection of a 3D nifti image. For each slice along
    the specified axis, the maximum or minimum intensity is computed from a number of
    neighboring slices around the current slice. The number of neighboring slices in
    each direction is defined by the size parameter.

    Parameters
    ----------
    file_in : str
        File name of input image.
    file_out : str
        File name of output image.
    size : int
        Number of included slices in each direction.
    axis : int, optional
        Axis along which to compute the intensity projection.
    mode : str, optional
        Mode (min, max).

    Returns
    -------
    None.

    """
    if axis < 0 or axis > 2:
        raise ValueError("Axis must be between 0 and 2.")

    # load input image
    data = nb.load(file_in)
    dim = data.header["dim"][axis + 1]
    arr = data.get_fdata()

    arr = np.moveaxis(arr, axis, 0)
    res = np.zeros_like(arr)
    for i in range(dim):
        index_low = i - size
        if index_low < 0:
            index_low = 0

        index_high = i + size
        if index_high > dim:
            index_high = dim

        if mode == "min":
            res[i, :, :] = np.min(arr[index_low:index_high, :, :], axis=0)
        elif mode == "max":
            res[i, :, :] = np.max(arr[index_low:index_high, :, :], axis=0)
        else:
            raise ValueError("Mode not supported")

    res = np.moveaxis(res, 0, axis)
    output = nb.Nifti1Image(res, data.affine, data.header)
    nb.save(output, file_out)


def remove_nans(file_in, file_out):
    """Reads a nifti volume and sets all nans to zero.

    Parameters
    ----------
    file_in : str
        Filename of input volume.
    file_out : str
        Filename of output volume.

    Returns
    -------
    None.

    """
    # load data
    data_img = nb.load(file_in)
    data_array = data_img.get_fdata()

    # set nans to zero
    data_array[np.isnan(data_array)] = 0

    # write output data
    output = nb.Nifti1Image(data_array, data_img.affine, data_img.header)
    nb.save(output, file_out)


def multiply_images(file1, file2, file_out):
    """This script does voxewise multiplication of two images. The output image takes
    the header information of the first image.

    Parameters
    ----------
    file1 : str
        Filename of first input image.
    file2 : str
        Filename of second input image.
    file_out : str
        Filename of the output image.

    Returns
    -------
    None.

    """
    # load both images
    file1_img = nb.load(file1)
    file2_img = nb.load(file2)
    file1_array = file1_img.get_fdata()
    file2_array = file2_img.get_fdata()

    # multiply both images
    file_out_array = file1_array * file2_array

    # write output image
    output = nb.Nifti1Image(file_out_array, file1_img.affine, file1_img.header)
    nb.save(output, file_out)


def average_layer(img_input, path_output, basename_output, mode="mean"):
    """This averages data across different layers or sessions. Input arrays should be in
    mgh format. The output gets the suffix of the chosen mode.

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
    for i, img in enumerate(img_input):
        data_res[:, i] = nb.load(img).get_fdata()[:, 0, 0]

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
    nb.save(output, os.path.join(path_output, basename_output + "_" + mode + ".mgh"))
