# -*- coding: utf-8 -*-

# external inputs
import nibabel as nb


def multiply_images(file1, file2, fileOUT):
    """ Multiply images

    This script does voxewise multiplication of two images. The output image 
    takes the header information of the first image.    

    Parameters
    ----------
    file1 : str
        Filename of first input image.
    file2 : str
        Filename of second input image.
    fileOUT : str
        Filename of the output image.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 01-11-2018             
    Last modified: 12-10-2020
    
    """
    
    # load both images
    file1_img = nb.load(file1)
    file2_img = nb.load(file2)
    file1_array = file1_img.get_fdata()
    file2_array = file2_img.get_fdata()
    
    # multiply both images
    fileOUT_array = file1_array * file2_array
    
    # write output image
    output = nb.Nifti1Image(fileOUT_array, file1_img.affine, file1_img.header)
    nb.save(output,fileOUT)
