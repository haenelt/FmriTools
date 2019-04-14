def get_acorr(input, write_output=False, path_output="", name_output=""):
    """
    This function computes a normalized autocorrelation of a 2D numpy array. The result is saved as 
    nifti image. The use of the scipy fftconvolve function is inspired by https://stackoverflow.com/
    questions/1100100/fft-based-2d-convolution-and-correlation-in-python.
    Inputs:
        *input: 2D nifti input array.
        *write_output: if output is written as nifti file (boolean).
        *path_output: path where output is saved.
        *name_output: basename of output image.
    Outputs:
        *array_corr: normalized autocorrelation array.
        
    created by Daniel Haenelt
    Date created: 11.04.2019
    Last modified: 12.04.2019
    """
    import os
    import numpy as np
    import nibabel as nb
    from scipy.signal import fftconvolve
    
    # normalize input array
    array1 = ( input - np.mean(input) ) / ( np.std(input) *np.shape(input)[0]*np.shape(input)[1] ) 
    array2 = ( input - np.mean(input) ) / ( np.std(input) ) 

    # compute autocorrelation   
    array_corr = fftconvolve(array1, array2[::-1, ::-1], mode="same")
    
    # normalize output
    array_corr = array_corr / np.max(array_corr)
    
    # write nifti
    if write_output is True:
        output = nb.Nifti1Image(array_corr, np.eye(4), nb.Nifti1Header())
        nb.save(output, os.path.join(path_output,name_output + "_nac.nii"))
    
    return array_corr
