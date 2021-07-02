# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb
from scipy.signal import fftconvolve


def get_acorr(arr, write_output=False, path_output="", name_output=""):
    """Get acorr.

    This function computes a normalized autocorrelation of a 2D numpy array. The 
    result is optionally saved as nifti image. The use of the scipy fftconvolve 
    function is inspired by [1]. The output autocorrelation is normalized to the 
    interval [0,1].    

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
