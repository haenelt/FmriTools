# -*- coding: utf-8 -*-

import nibabel as nb


def multiply_images(file1, file2, file_out):
    """Multiply images.

    This script does voxewise multiplication of two images. The output image
    takes the header information of the first image.

    Parameters
    ----------
    file1 : str
        Filename of first input image.
    file2 : str
        Filename of second input image.
    file_out : str
        Filename of the output image.

    Returns
    -------
    None.

    """

    # load both images
    file1_img = nb.load(file1)
    file2_img = nb.load(file2)
    file1_array = file1_img.get_fdata()
    file2_array = file2_img.get_fdata()

    # multiply both images
    file_out_array = file1_array * file2_array

    # write output image
    output = nb.Nifti1Image(file_out_array, file1_img.affine, file1_img.header)
    nb.save(output, file_out)
