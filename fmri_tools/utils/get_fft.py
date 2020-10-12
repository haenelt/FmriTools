# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb
from numpy.fft import fft2, fftshift
from numpy.random import shuffle

    
def get_fft(input, write_output=False, path_output="", name_output="", 
            normalization=False, N=1000):
    """
    This function computes the power spectrum of a 2D numpy array. The result is 
    saved as nifti image. Optionally, the FFT spectrum can be normalized. To do 
    so, the input is shuffled N times and the 50th percentile of the power 
    spectrum of the shuffled input array is computed. The mean 50th percentile 
    of all shuffles is subtracted from the initial power spectrum and negative 
    values are set to zero. Additionally, the resulting spectrum is normalized 
    by its maximum power.
    Inputs:
        *input: 2D nifti input array.
        *write_output: if output is written as nifti file (boolean).
        *path_output: path where output is saved.
        *name_output: basename of output image.
        *normalization: indicate if power spectrum is normalized (boolean).
        *N: number of shuffled for normalization.
    Outputs:
        *array_fft: fourier transformed array.
        
    created by Daniel Haenelt
    Date created: 11-04-2019
    Last modified: 12-10-2020
    """
    
    # copy input data
    data = input.copy()
    
    # compute autocorrelation
    array_fft = np.abs(fftshift(fft2(data)))

    shuffle_mean = []
    if normalization is True:
        for i in range(N):
            # shuffle input array
            shuffle(data) # shuffle along first axis
            data = data.T # transpose array
            shuffle(data) # shuffle along second axis
            data = data.T            
        
            # get fft
            array_shuffle_fft = np.abs(fftshift(fft2(data)))
            
            # get 50th percentile
            shuffle_mean.append(np.percentile(array_shuffle_fft,50))
        
        # get mean of percentiles
        array_50 = np.mean(shuffle_mean)
        
        # subtract 50th percentile from fft and divide by maximum power
        array_fft = array_fft - array_50
        array_fft[array_fft < 0] = 0
        array_fft = array_fft / np.max(array_fft)
        
    # write nifti
    if write_output is True:
        output = nb.Nifti1Image(array_fft, np.eye(4), nb.Nifti1Header())
        nb.save(output, os.path.join(path_output,name_output + "_fft.nii"))
        
    return array_fft
