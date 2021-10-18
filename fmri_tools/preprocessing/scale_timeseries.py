# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb

# local inputs
from ..io.get_filename import get_filename


def scale_timeseries(file_in, cutoff=50, prefix="p"):
    """Scale time series.

    The function converts a time series to percent signal change. Each voxel is
    divided by its temporal mean and multiplied with 100.

    Parameters
    ----------
    file_in : str
        Filename of nifti time series.
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

    # center around zero
    arr -= 100

    # write output
    output = nb.Nifti1Image(arr, data.affine, data.header)
    nb.save(output, os.path.join(path_file, prefix + name_file + ext_file))
