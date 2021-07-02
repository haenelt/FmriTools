# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb

# local inputs
from ..io.get_filename import get_filename


def get_tsnr(file_in, tsnr_max=200, write_output=False, path_output=""):
    """Get tSNR.
    
    This function computes the tsnr of one time series.

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
        nb.save(data_img, os.path.join(path_output, "tsnr_"+file+ext))
    
    return data_tsnr_array
