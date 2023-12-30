# -*- coding: utf-8 -*-

import os

import nibabel as nb
import numpy as np


def get_gaussian(input_file, path_output, mu=0, sigma=1):
    """Get gaussian.

    Creates a NIfTI volume with Gaussian noise as voxel intensities and saves a
    nifti volume with matrix size and header information taken from the input
    image.

    Parameters
    ----------
    input_file : str
        Input image.
    path_output : str
        Path where output file is written.
    mu : float, optional
        Mean of Gaussian noise. The default is 0.
    sigma : float, optional
        Standard deviation of Gaussian noise. The default is 1.

    Returns
    -------
    None.

    """

    # load img
    img = nb.load(input_file)

    # new image size
    x_size = img.shape[0]
    y_size = img.shape[1]
    z_size = img.shape[2]

    # random Gaussian noise for each pixel
    img_data = np.random.normal(mu, sigma, x_size * y_size * z_size)
    img_data = np.reshape(img_data, (x_size, y_size, z_size))

    # write output image
    newimg = nb.Nifti1Image(img_data, img.affine, img.header)
    nb.save(newimg, os.path.join(path_output, "gaussian_noise.nii"))
