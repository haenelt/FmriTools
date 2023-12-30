# -*- coding: utf-8 -*-

import os

import nibabel as nb
import numpy as np


def get_surface_voxel(file_in, path_output):
    """Get surface voxel.

    Compute surface voxels for to illustrate curvature dependence in surface
    representations as introduced in [1].

    Parameters
    ----------
    file_in : str
        Input image.
    path_output : str
        Path where to write output image.

    Returns
    -------
    None.

    References
    -------
    .. [1] Kay, K, et al. A critical assessment of data quality and venous
    effects in ultra-high-resolution fMRI, bioRxiv, 1--45 (2018).

    """

    # parameters
    x_calc = 1
    y_calc = 1
    z_calc = 1
    x_max = 1
    y_max = 2
    z_max = 4

    # make subfolders
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # load img
    img = nb.load(file_in)

    # img range
    x_size = img.shape[0]
    y_size = img.shape[1]
    z_size = img.shape[2]

    # compute surface voxels
    img_res = np.zeros([x_size, y_size, z_size])
    for i in range(x_size):
        for j in range(y_size):
            for k in range(z_size):
                if np.mod(i, 2) != 0 and x_calc != 0:
                    img_res[i, j, k] += x_max
                if np.mod(j, 2) != 0 and y_calc != 0:
                    img_res[i, j, k] += y_max
                if np.mod(k, 2) != 0 and z_calc != 0:
                    img_res[i, j, k] += z_max

    # write output image
    newimg = nb.Nifti1Image(img_res, img.affine, img.header)
    nb.save(newimg, os.path.join(path_output, "surface_voxel.nii"))
