# -*- coding: utf-8 -*-
"""Different image manipulation tools."""

import os
import subprocess
import sys

import nibabel as nb
import numpy as np

from ..io.filename import get_filename
from ..utils.bias import remove_bias_ants

__all__ = [
    "estimate_t1w",
    "mip",
    "remove_nans",
    "multiply_images",
    "average_layer",
    "laminar_profile",
    "volume_threshold",
    "mean_image",
    "create_series",
]


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


def laminar_profile(file_in, path_output, hemi, name_output, mode):
    """This function computes the laminar profile (mean, median, max, min) of data
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
    nb.save(output, os.path.join(path_output, f"{hemi}.{name_output}.mgh"))


def volume_threshold(file_in, prefix, data_max):
    """Takes a nifti volume and sets a maximum threshold value. All values above the
    threshold are replaced by the threshold value. The output image is saved in the same
    folder. A prefix is added to the file name.

    Parameters
    ----------
    file_in : str
        Path of input image.
    prefix : str
        Defined prefix for the output image.
    data_max : float
        Set maximum threshold value.

    Returns
    -------
    None.

    """
    # load data
    data_img = nb.load(file_in)
    data_array = data_img.get_fdata()

    # set maximum threshold
    data_array[data_array > data_max] = data_max

    # write output data
    filename_out = os.path.join(
        os.path.dirname(file_in), prefix + os.path.basename(file_in)
    )
    output = nb.Nifti1Image(data_array, data_img.affine, data_img.header)
    nb.save(output, filename_out)


def mean_image(file_in, file_out):
    """Compute temporal mean of time series using FSL.

    Parameters
    ----------
    file_in : str
        File name of input time series.
    file_out : str
        File name of temporal mean image.
    """
    # get filename
    path_out, _, _ = get_filename(file_out)

    # make output folder
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    command = "fslmaths"
    command += f" {file_in}"
    command += " -Tmean"
    command += f" {file_out}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")


def create_series(file_in, path_out, name_output):
    """This function creates a 4D nifti time series from a set of 3D nifti files.

    Parameters
    ----------
    file_in : list
        Array of filenames containing single 3D nifti volumes.
    path_out : str
        Path where output is saved.
    name_output : str
        Basename of output 4D nifti file.

    Returns
    -------
    None.

    """
    # load 3D nifti to get array size
    data = nb.load(os.path.join(file_in[0]))
    data.header["dim"][0] = 4
    data.header["dim"][4] = len(file_in)

    res = np.zeros(data.header["dim"][1:5])
    for i in range(len(file_in)):
        img = nb.load(os.path.join(file_in[i])).get_fdata()
        res[:, :, :, i] = img

    output = nb.Nifti1Image(res, data.affine, data.header)
    nb.save(output, os.path.join(path_out, name_output + ".nii"))
