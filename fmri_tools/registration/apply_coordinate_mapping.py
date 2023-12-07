# -*- coding: utf-8 -*-

import os

import nibabel as nb
import numpy as np
from fmri_tools.utils.interpolation import linear_interpolation3d, nn_interpolation3d

# linear and nearest neighbor sampling functions
_sampler = {"linear": linear_interpolation3d, "nearest": nn_interpolation3d}


def _set_min(arr, min_val):
    """Remove coordinates below the matrix size."""
    arr[np.floor(arr) < min_val] = None
    return arr


def _set_max(arr, max_val):
    """Remove coordinates above the matrix size."""
    arr[np.ceil(arr) > max_val] = None
    return arr


def apply_coordinate_mapping(file_in, cmap_in, file_out, interpolation="linear"):
    """Apply coordinate mapping.

    This function applies a coordinate mapping to a volume.

    Parameters
    ----------
    file_in : str
        Filename of input volume.
    cmap_in : str
        Filename of coordinate mapping.
    file_out : str
        Filename of output volume.
    interpolation : str, optional
        Interpolation type (linear or nearest). The default is "linear".

    Returns
    -------
    niimg
        Transformed volume.

    """

    # make output folder
    path_output = os.path.dirname(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # load data
    data = nb.load(file_in)
    arr = data.get_fdata()
    cmap = nb.load(cmap_in)
    arr_c = cmap.get_fdata()

    # get source and target image dimensions
    x_dim_source = data.header["dim"][1]
    y_dim_source = data.header["dim"][2]
    z_dim_source = data.header["dim"][3]
    x_dim_target = cmap.header["dim"][1]
    y_dim_target = cmap.header["dim"][2]
    z_dim_target = cmap.header["dim"][3]

    # get mapping coordinates
    arr_c_x = arr_c[:, :, :, 0]
    arr_c_y = arr_c[:, :, :, 1]
    arr_c_z = arr_c[:, :, :, 2]

    # flatten
    arr_c_x = arr_c_x.flatten()
    arr_c_y = arr_c_y.flatten()
    arr_c_z = arr_c_z.flatten()

    # remove boundary coordinates
    arr_c_x = _set_min(arr_c_x, min_val=0.0)
    arr_c_y = _set_min(arr_c_y, min_val=0.0)
    arr_c_z = _set_min(arr_c_z, min_val=0.0)
    arr_c_x = _set_max(arr_c_x, max_val=x_dim_source - 1)
    arr_c_y = _set_max(arr_c_y, max_val=y_dim_source - 1)
    arr_c_z = _set_max(arr_c_z, max_val=z_dim_source - 1)

    # get coordinates to keep
    arr_sum = arr_c_x + arr_c_y + arr_c_z
    ind_ignore = np.where(np.isnan(arr_sum))[0]
    ind_keep = np.arange(len(arr_c_x))
    ind_keep = np.array(list(set(ind_keep).difference(set(ind_ignore))))

    # only use kept coordinates for interpolation
    arr_c_x = arr_c_x[ind_keep]
    arr_c_y = arr_c_y[ind_keep]
    arr_c_z = arr_c_z[ind_keep]

    # do the interpolation
    arr_sampled = _sampler[interpolation](arr_c_x, arr_c_y, arr_c_z, arr)

    # reshape to output array
    res = np.zeros_like(arr_sum)
    res[ind_keep] = arr_sampled
    res = np.reshape(res, (x_dim_target, y_dim_target, z_dim_target))

    output = nb.Nifti1Image(res, cmap.affine, cmap.header)
    nb.save(output, file_out)

    return output
