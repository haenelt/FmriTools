# -*- coding: utf-8 -*-

import nibabel as nb
import numpy as np


def calc_mip(file_in, file_out, size, axis=2, mode="min"):
    """Minimum/Maximum intensity projection of a 3D nifti image. For each slice along
    the specified axis, the maximum or minimum intensity is computed from a number of
    neighboring slices around the current slice. The number of neighboring slices in
    each direction is defined by the size parameter.

    Parameters
    ----------
    file_in : str
        File name of input image.
    file_out : str
        File name of output image.
    size : int
        Number of included slices in each direction.
    axis : int, optional
        Axis along which to compute the intensity projection.
    mode : str, optional
        Mode (min, max).

    Returns
    -------
    None.

    """

    if axis < 0 or axis > 2:
        raise ValueError("Axis must be between 0 and 2.")

    # load input image
    data = nb.load(file_in)
    dim = data.header["dim"][axis + 1]
    arr = data.get_fdata()

    arr = np.moveaxis(arr, axis, 0)
    res = np.zeros_like(arr)
    for i in range(dim):
        index_low = i - size
        if index_low < 0:
            index_low = 0

        index_high = i + size
        if index_high > dim:
            index_high = dim

        if mode == "min":
            res[i, :, :] = np.min(arr[index_low:index_high, :, :], axis=0)
        elif mode == "max":
            res[i, :, :] = np.max(arr[index_low:index_high, :, :], axis=0)
        else:
            raise ValueError("Mode not supported")

    res = np.moveaxis(res, 0, axis)
    output = nb.Nifti1Image(res, data.affine, data.header)
    nb.save(output, file_out)
