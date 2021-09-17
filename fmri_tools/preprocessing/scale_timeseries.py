# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb

# local inputs
from ..io.get_filename import get_filename


def scale_timeseries(file_in, prefix="p"):
    """Scale time series.

    The function converts a time series to percent signal change. Each voxel is
    divided by its temporal mean and multiplied with 100.

    Parameters
    ----------
    file_in : str
        Filename of nifti time series.
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

    # load array with appended volumes
    data_array = data.get_fdata()
    mean_array = np.mean(data_array, axis=3)

    # compute percent signal change
    data_array = data_array / mean_array * 100

    # write output
    output = nb.Nifti1Image(data_array, data.affine, data.header)
    output.set_data_dtype(np.float)
    nb.save(output, os.path.join(path_file, prefix + name_file + ext_file))
