# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb

# local inputs
from ..io.get_filename import get_filename

__all__ = ['scale_timeseries', 'demean_timeseries']


def scale_timeseries(file_in, cutoff=50, prefix="p"):
    """Scale time series.

    The function converts a time series to percent signal change. Each voxel is
    divided by its temporal mean and multiplied with 100.

    Parameters
    ----------
    file_in : str
        Filename.
    cutoff : float, optional
        Threshold data which exceed the cutoff percentage.
    prefix : str, optional
        Prefix of output time series basename. The default is "p".

    Returns
    -------
    None.
    
    """

    # get filename
    path_file, name_file, ext_file = get_filename(file_in)

    # load data
    data = nb.load(file_in)
    nt = data.header["dim"][4]

    # load array with appended volumes
    arr = data.get_fdata()
    arr_mean = np.mean(arr, axis=3)
    arr_mean = np.repeat(arr_mean[:, :, :, np.newaxis], nt, axis=3)

    # compute percent signal change
    arr = arr / arr_mean * 100
    
    # threshold data
    if cutoff:
        arr_max = np.max(arr, axis=3)
        arr_min = np.min(arr, axis=3)
    
        arr[arr_max > 100 + cutoff, :] = 100 + cutoff
        arr[arr_min < 100 - cutoff, :] = 100 - cutoff

    # write output
    output = nb.Nifti1Image(arr, data.affine, data.header)
    nb.save(output, os.path.join(path_file, prefix + name_file + ext_file))


def demean_timeseries(img_input, path_output="", name_output="",
                      write_output=False):
    """De-mean time series.

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
        data_array[:, :, :, i] = (data_array[:, :, :,
                                  i] - data_mean) / data_mean * 100

    # write output
    output = nb.Nifti1Image(data_array, img_input.affine, img_input.header)
    if write_output:
        nb.save(output,
                os.path.join(path_output, "demean_" + name_output + ".nii"))

    return output
