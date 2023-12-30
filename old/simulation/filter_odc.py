# -*- coding: utf-8 -*-

import numpy as np
from numpy.fft import fftshift


def filter_odc_2d(nx, ny, fovx, fovy, rho, delta, epsilon, theta):
    """Filter ODC 2D.

    This function defines the 2D ocular dominance column filter (anisotropic
    band-pass filter) in spatial frequency space which is taken from [1]. The
    filter is normalised in order to have the same variance as the white noise
    source image.

    Parameters
    ----------
    nx : int
        Array size in x-direction.
    ny : int
        Array size in y-direction.
    fovx : float
        Field of view in x-direction (mm).
    fovy : float
        Field of view in y-direction (mm).
    rho : float
        Main spatial frequency determining columnar width in cycles/mm.
    delta : float
        Variations orthogonal to ODC bands (irregularity).
    epsilon : float
        Variations parallel to ODC bands (branchiness).
    theta : float
        Rotation angle in deg.

    Returns
    -------
    f : ndarray
        Band-pass filter array in spatial frequency space.

    References
    -------
    .. [1] Chaimow, D, et al. Spatial specificity of the functional MRI blood
    oxygenation response relative to neuronal activity, Neuroimage 164, 32--47
    (2018).

    """

    # get maximum k-space coordinate in x- and y-direction
    kx_max = nx / (2 * fovx)
    ky_max = ny / (2 * fovy)

    # get k-space axes
    kx = np.linspace(-kx_max, kx_max, nx)
    ky = np.linspace(-ky_max, ky_max, ny)

    # define two-dimensional k-space grid
    kx_grid, ky_grid = np.meshgrid(kx, ky)

    # convert theta from deg to rad
    theta = np.radians(theta)

    # convert to polar coordinates
    k_r = np.sqrt(kx_grid**2 + ky_grid**2)
    k_pol = np.arctan2(ky_grid, kx_grid)

    # ODC filter
    f_rad1 = np.exp(-((k_r - rho) ** 2) / (2 * delta**2))
    f_rad2 = np.exp(-((k_r + rho) ** 2) / (2 * delta**2))

    f_ang1 = np.exp(np.cos(k_pol - theta) / epsilon**2)
    f_ang2 = np.exp(-np.cos(k_pol - theta) / epsilon**2)

    f = (f_rad1 + f_rad2) * (f_ang1 + f_ang2)

    # normalise filter
    c = np.sqrt(np.sum(f**2) / np.size(f))

    f = f / c

    # shift zero k-space line to the left border to match the numpy fft
    # convention
    f = fftshift(f)

    return f


def filter_odc_1d(n, fov, rho, delta):
    """Filter ODC 1D.

    This function defines the 1D ocular dominance column filter (anisotropic
    band-pass filter) in spatial frequency space which is taken from [1]. The
    filter is normalised in order to have the same variance as the white noise
    source image.

    Parameters
    ----------
    n : int
        Array size.
    fov : float
        Field of view (mm).
    rho : float
        Main spatial frequency determining columnar width in cycles/mm.
    delta : float
        Variations orthogonal to ODC bands (irregularity).

    Returns
    -------
    f : ndarray
        Band-pass filter array in spatial frequency space.

    References
    -------
    .. [1] Chaimow, D, et al. Spatial specificity of the functional MRI blood
    oxygenation response relative to neuronal activity, Neuroimage 164, 32--47
    (2018).

    """

    # get maximum k-space coordinate in x- and y-direction
    k_max = n / (2 * fov)

    # get k-space axes
    k = np.linspace(-k_max, k_max, n)

    # ODC filter
    f_rad1 = np.exp(-((k - rho) ** 2) / (2 * delta**2))
    f_rad2 = np.exp(-((k + rho) ** 2) / (2 * delta**2))

    f = f_rad1 + f_rad2

    # normalise filter
    c = np.sqrt(np.sum(f**2) / np.size(f))

    f = f / c

    # shift zero k-space line to the left border to match the numpy fft
    # convention
    f = fftshift(f)

    return f
