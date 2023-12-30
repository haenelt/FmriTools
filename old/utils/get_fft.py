# -*- coding: utf-8 -*-

import os

import nibabel as nb
import numpy as np
from numpy.fft import fft2, fftshift
from numpy.random import shuffle


def get_fft(
    arr, write_output=False, path_output="", name_output="", normalization=False, n=1000
):
    """Get FFT.

    This function computes the power spectrum of a 2D numpy array. The result is
    saved as nifti image. Optionally, the FFT spectrum can be normalized. To do
    so, the input is shuffled N times and the 50th percentile of the power
    spectrum of the shuffled input array is computed. The mean 50th percentile
    of all shuffles is subtracted from the initial power spectrum and negative
    values are set to zero. Additionally, the resulting spectrum is normalized
    by its maximum power.

    Parameters
    ----------
    arr : ndarray
        2D nifti input array.
    write_output : bool, optional
        If output is written as nifti file. The default is False.
    path_output : str, optional
        Path where output is saved. The default is "".
    name_output : str, optional
        Basename of output image. The default is "".
    normalization : bool, optional
        Indicate if power spectrum is normalized. The default is False.
    n : int, optional
        Number of shuffled for normalization. The default is 1000.

    Returns
    -------
    array_fft : ndarray
        Fourier transformed array.

    """

    # copy input data
    data = arr.copy()

    # compute autocorrelation
    array_fft = np.abs(fftshift(fft2(data)))

    shuffle_mean = []
    if normalization is True:
        for i in range(n):
            # shuffle input array
            shuffle(data)  # shuffle along first axis
            data = data.T  # transpose array
            shuffle(data)  # shuffle along second axis
            data = data.T

            # get fft
            array_shuffle_fft = np.abs(fftshift(fft2(data)))

            # get 50th percentile
            shuffle_mean.append(np.percentile(array_shuffle_fft, 50))

        # get mean of percentiles
        array_50 = np.mean(shuffle_mean)

        # subtract 50th percentile from fft and divide by maximum power
        array_fft = array_fft - array_50
        array_fft[array_fft < 0] = 0
        array_fft = array_fft / np.max(array_fft)

    # write nifti
    if write_output is True:
        output = nb.Nifti1Image(array_fft, np.eye(4), nb.Nifti1Header())
        nb.save(output, os.path.join(path_output, name_output + "_fft.nii"))

    return array_fft
