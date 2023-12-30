# -*- coding: utf-8 -*-

import os

import nibabel as nb


def volume_threshold(file_in, prefix, data_max):
    """Volume threshold.

    Takes a nifti volume and sets a maximum threshold value. All values above
    the threshold are replaced by the threshold value. The output image is saved
    in the same folder. A prefix is added to the file name.

    Parameters
    ----------
    file_in : str
        Path of input image.
    prefix : str
        Defined prefix for the output image.
    data_max : float
        Set maximum threshold value.

    Returns
    -------
    None.

    """

    # load data
    data_img = nb.load(file_in)
    data_array = data_img.get_fdata()

    # set maximum threshold
    data_array[data_array > data_max] = data_max

    # write output data
    filename_out = os.path.join(
        os.path.dirname(file_in), prefix + os.path.basename(file_in)
    )
    output = nb.Nifti1Image(data_array, data_img.affine, data_img.header)
    nb.save(output, filename_out)
