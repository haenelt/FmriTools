# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys

# external inputs
import numpy as np
import nibabel as nb


def get_mip(file_in, slice_dir, slice_start, slice_end, path_output,
            name_output):
    """Get MIP.
    
    Create minimum intensity projection between between slices of a nifti 
    volume.

    Parameters
    ----------
    file_in : str
        Input image.
    slice_dir : int
        Projection direction (0,1,2).
    slice_start : int
        Starting slice index for minimum intensity projection.
    slice_end : int
        Ending slice index for minimum intensity projection.
    path_output : str
        Path where to write output image.
    name_output : str
        Basename of output file.

    Raises
    ------
    ValueError
        If `slice_dir` is invalid.

    Returns
    -------
    None.

    """

    # make subfolders
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # load data
    data_img = nb.load(file_in)
    data_array = data_img.get_data()

    # new image size
    if slice_dir == 0:
        nx = data_img.shape[1]
        ny = data_img.shape[2]
    elif slice_dir == 1:
        nx = data_img.shape[0]
        ny = data_img.shape[2]
    elif slice_dir == 2:
        nx = data_img.shape[0]
        ny = data_img.shape[1]
    else:
        sys.exit("Choose valid slice direction!")

    # get minimum projection
    if slice_dir == 0:
        slice_array = data_array[slice_start:slice_end, :, :]
    elif slice_dir == 1:
        slice_array = data_array[:, slice_start:slice_end, :]
    elif slice_dir == 2:
        slice_array = data_array[:, :, slice_start:slice_end]
    else:
        raise ValueError("Unknown slice direction.")

    # mip
    mip = np.min(slice_array, axis=slice_dir)

    # output image matrix
    newimg = nb.Nifti1Image(mip, data_img.affine, data_img.header)
    newimg.header["dim"][0] = 2
    newimg.header["dim"][1] = nx
    newimg.header["dim"][2] = ny
    newimg.header["dim"][3] = 1
    newimg.header["dim"][4] = 1
    nb.save(newimg, os.path.join(path_output, name_output + "_mip.nii"))
