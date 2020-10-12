# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb
    
    
def mask_ana(t1, mask, background_bright=False):
    """
    This function masked an image with a corresponding binary mask by 
    multiplication. The masked volume is saved in the same folder as the input 
    image with the prefix p.
    Inputs:
        *t1: input anatomy.
        *mask: corresponding binary mask.
        *background_bright: set values outside mask to maximum value (boolean).
        
    created by Daniel Haenelt
    Date created: 13-02-2019   
    Last modified: 12-10-2020
    """

    # get path and filename of anatomy
    path_t1 = os.path.dirname(t1)
    if os.path.splitext(os.path.basename(t1))[1] == '.gz':
        name_t1 = os.path.splitext(os.path.splitext(os.path.basename(t1))[0])[0]
    else:
        name_t1 = os.path.splitext(os.path.basename(t1))[0]

    # load anatomy
    ana_img = nb.load(t1)
    ana_array = ana_img.get_fdata()
    
    # load mask
    mask_img = nb.load(mask)
    mask_array = mask_img.get_fdata()
    
    # multiply images
    masked_ana_array = ana_array * mask_array
    
    # set all outside mask values to maximum to mimic bright CSF values
    if background_bright:
        masked_ana_array[mask_array == 0] = np.max(ana_array)
    
    # write masked anatomy
    out_img = nb.Nifti1Image(masked_ana_array, ana_img.affine, ana_img.header)
    nb.save(out_img,os.path.join(path_t1,"p"+name_t1+".nii"))
    