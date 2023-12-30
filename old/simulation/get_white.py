# -*- coding: utf-8 -*-

import numpy as np


def get_white_2d(nx, ny, mu, sigma):
    """Get white 2D.

    Creates a 2D array with Gaussian noise.

    Parameters
    ----------
    nx : int
        Array size in x-direction.
    ny : int
        Array size in x-direction.
    mu : float
        Mean of the Gaussian distribution.
    sigma : float
        Standard deviation of the Gaussian distribution.

    Returns
    -------
    img : ndarray
        White noise array.

    """

    img = np.random.normal(mu, sigma, nx * ny)
    img = np.reshape(img, (nx, ny))

    return img


def get_white_1d(n, mu, sigma):
    """Get white 1D.

    Creates a 1D array with Gaussian noise.

    Parameters
    ----------
    n : int
        Array size.
    mu : float
        Mean of the Gaussian distribution.
    sigma : float
        Standard deviation of the Gaussian distribution.

    Returns
    -------
    img : ndarray
        White noise array.

    """

    img = np.random.normal(mu, sigma, n)

    return img
