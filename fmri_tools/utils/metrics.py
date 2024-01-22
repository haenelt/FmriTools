# -*- coding: utf-8 -*-

import copy
import os

import nibabel as nb
import numpy as np
from scipy.signal import fftconvolve

from ..io.filename import get_filename

__all__ = ["calc_tsnr", "calc_mean", "calc_mean4d", "calc_std", "calc_acorr"]


def calc_tsnr(file_in, tsnr_max=200, write_output=False, path_output=""):
    """This function computes the tsnr of one time series.

    Parameters
    ----------
    file_in : str
        Input time series.
    tsnr_max : TYPE, optional
        Threshold unrealistic high tsnr values (applied if set > 0). The default
        is 200.
    write_output : bool, optional
        Write output nifti file. The default is False.
    path_output : str, optional
        Path where to save mean image The default is "".

    Returns
    -------
    data_tsnr_array : ndarray
        TSNR array.

    """
    # make subfolders
    if write_output and not os.path.exists(path_output):
        os.makedirs(path_output)

    # get filename
    _, file, ext = get_filename(file_in)

    # load time series
    data_img = nb.load(file_in)
    data_array = data_img.get_fdata()

    # get mean and std
    data_mean_array = np.mean(data_array, axis=3)
    data_std_array = np.std(data_array, axis=3)
    data_std_array[data_std_array == 0] = np.nan  # set zeroes to nan

    # get tsnr of time series
    data_tsnr_array = data_mean_array / data_std_array
    data_tsnr_array[np.isnan(data_tsnr_array)] = 0

    # threshold tsnr
    if tsnr_max:
        data_tsnr_array[data_tsnr_array > tsnr_max] = tsnr_max

    # write output
    if write_output:
        data_img.header["dim"][0] = 3
        data_img.header["dim"][4] = 1

        data_img = nb.Nifti1Image(data_tsnr_array, data_img.affine, data_img.header)
        nb.save(data_img, os.path.join(path_output, "tsnr_" + file + ext))

    return data_tsnr_array


def calc_mean(file_in, path_output, name_output, method="mean"):
    """This function computes the mean image of one or more time series.

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
        for i, _ in enumerate(file_in):
            data_img = nb.load(file_in[i])
            t_size += data_img.header["dim"][4]

        data_all_array = np.zeros([x_size, y_size, z_size, t_size])
        t_start = 0
        for i, _ in enumerate(file_in):
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


def calc_mean4d(file_in, path_output="", name_output="", write_output=False):
    """This function computes the mean time series of one or more time series.

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
    for i, _ in enumerate(file_in):
        res_array += nb.load(file_in[i]).get_fdata()

    res_array = res_array / len(file_in)

    # write mean time series
    output = nb.Nifti1Image(res_array, data_img.affine, data_img.header)
    if write_output:
        nb.save(output, os.path.join(path_output, "mean_" + name_output + ".nii"))

    return output


def calc_std(file_in, path_output, name_output, set_outlier=None):
    """Get std.

    This function computes the standard deviation of one or more time series.

    Parameters
    ----------
    file_in : str
        Single file or list of files.
    path_output : str
        Path where to save mean image
    name_output : str
        Output file name without file extension.
    set_outlier : float, optional
        Can be nan, zero or None. The default is None.
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
        for i, _ in enumerate(file_in):
            data_img = nb.load(file_in[i])
            t_size += data_img.header["dim"][4]

        data_all_array = np.zeros([x_size, y_size, z_size, t_size])
        t_start = 0
        for i, _ in enumerate(file_in):
            data_img = nb.load(file_in[i])
            data_array = data_img.get_fdata()
            t_end = data_img.header["dim"][4] + t_start
            data_all_array[:, :, :, t_start:t_end] = data_array
            t_start = copy.deepcopy(t_end)

    else:
        # load data
        data_img = nb.load(file_in)
        data_all_array = data_img.get_fdata()

    # calculate std
    data_std_array = np.std(data_all_array, axis=3)

    if set_outlier == "nan":
        data_std_array[data_std_array == 0] = np.nan  # set zeroes to nan
    elif set_outlier == "zero":
        data_std_array[data_std_array == 0] = 0

    # write output
    if len(np.shape(file_in)) > 0:
        data_img = nb.load(file_in[0])
    else:
        data_img = nb.load(file_in)
    data_img.header["dim"][0] = 3
    data_img.header["dim"][4] = 1

    # write mean image
    mean_img = nb.Nifti1Image(data_std_array, data_img.affine, data_img.header)
    nb.save(mean_img, os.path.join(path_output, "std_" + name_output + ".nii"))


def calc_acorr(arr, write_output=False, path_output="", name_output=""):
    """This function computes a normalized autocorrelation of a 2D numpy array. The
    result is optionally saved as nifti image. The use of the scipy fftconvolve function
    is inspired by [1]. The output autocorrelation is normalized to the interval [0,1].

    Parameters
    ----------
    arr : ndarray
        2D nifti input array.
    write_output : bool, optional
        If output is written as nifti fil. The default is False.
    path_output : str, optional
        Path where output is saved. The default is "".
    name_output : str, optional
        Basename of output image. The default is "".

    Returns
    -------
    array_corr : ndarray
        Normalized autocorrelation array.

    References
    -------
    .. [1] https://stackoverflow.com/questions/1100100/fft-based-2d-convolution-
    and-correlation-in-python

    """
    # normalize input array
    array1 = (arr - np.mean(arr)) / (np.std(arr) * np.shape(arr)[0] * np.shape(arr)[1])
    array2 = (arr - np.mean(arr)) / (np.std(arr))

    # compute autocorrelation
    array_corr = fftconvolve(array1, array2[::-1, ::-1], mode="same")

    # normalize output
    array_corr = array_corr / np.max(array_corr)

    # write nifti
    if write_output is True:
        output = nb.Nifti1Image(array_corr, np.eye(4), nb.Nifti1Header())
        nb.save(output, os.path.join(path_output, name_output + "_nac.nii"))

    return array_corr
