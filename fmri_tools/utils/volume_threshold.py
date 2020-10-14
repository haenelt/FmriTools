# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import nibabel as nb

    
def volume_threshold(filename, prefix, data_max):
    """ Volume threshold

    Takes a nifti volume and sets a maximum threshold value. All values above 
    the threshold are replaced by the threshold value. The output image is saved 
    in the same folder. A prefix is added to the file name.    

    Parameters
    ----------
    filename : str
        Path of input image.
    prefix : str
        Defined prefix for the output image.
    data_max : float
        Set maximum threshold value.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 01-11-2018             
    Last modified: 12-10-2020

    """
    
    # load data
    data_img = nb.load(filename)
    data_array = data_img.get_fdata()
    
    # set maximum threshold
    data_array[data_array > data_max] = data_max
    
    # write output data
    filenameOUT = os.path.join(os.path.dirname(filename),prefix+os.path.basename(filename))
    output = nb.Nifti1Image(data_array, data_img.affine, data_img.header)
    nb.save(output,filenameOUT)
    