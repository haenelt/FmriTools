# -*- coding: utf-8 -*-
"""Simulation code for the organization of cortical columns."""

import math
import sys

import nibabel as nb
import numpy as np
from numpy.fft import fft, fft2, fftshift, ifft, ifft2
from scipy.interpolate import griddata
from scipy.stats import pearsonr

__all__ = [
    "odc_2d",
    "odc_1d",
    "pattern_2d",
    "pattern_1d",
    "filter_bold_2d",
    "filter_bold_1d",
    "filter_odc_2d",
    "filter_odc_1d",
    "filter_sigmoid",
    "get_white_2d",
    "get_white_1d",
    "mask_pattern_2d",
    "mask_pattern_1d",
    "regrid_zero_2d",
    "regrid_zero_1d",
    "regrid_2d",
    "regrid_1d",
    "pattern_corr",
]


def odc_2d(
    nx_sim=1024,
    ny_sim=1024,
    fov_x=20,
    fov_y=20,
    nx_mri=100,
    ny_mri=100,
    rho=0.5,
    delta=0.3,
    epsilon=0.4,
    theta=0,
    alpha=4,
    beta=0.05,
    fwhm_bold=1.02,
    fwhm_noise=0.001,
    a_mask=1000,
    b_mask=1000,
    alpha_mask=0,
    path_white=None,
):
    """This function generates realistic ocular dominance patterns according to the
    model proposed by [1] in 2D. Implementation follows [2] and most of the default
    values for the columnar pattern are taken from this publication. First, a white
    noise pattern is defined and an anisotropic band-pass filter and a non-linear
    sigmoidal filter are applied. This neural map is converted to a BOLD response map by
    applying a Gaussian filter. White noise is added to simulate uncorrelated
    measurement noise to the signal. The MRI sampling procedure is considered by only
    inner k-space lines of the spatial frequency representation of the BOLD reseponse
    map.

    Parameters
    ----------
    nx_sim : int, optional
        Array size of the simulated patch in x-direction. The default is 1024.
    ny_sim : int, optional
        Array size of the simulated patch in y-direction. The default is 1024.
    fov_x : float, optional
        Field of view in x-direction (mm). The default is 20.
    fov_y : float, optional
        Field of view in y-direction (mm). The default is 20.
    nx_mri : int, optional
        Array size of the MR image in x-direction. The default is 100.
    ny_mri : int, optional
        Array size of the MR image in y-direction. The default is 100.
    rho : float, optional
        Main spatial frequency determining columnar width in cycles/mm. The
        default is 0.5.
    delta : float, optional
        Variations orthogonal to ODC bands in cycles/mm (irregularity). The
        default is 0.3.
    epsilon : float, optional
        Variations parallel to ODC bands in cycles/mm (branchiness). The default
        is 0.4.
    theta : float, optional
        Orientation of the columnar pattern in deg. The default is 0.
    alpha : float, optional
        Sharpness parameter of the sigmoidal filter. The default is 4.
    beta : float, optional
        Maximal BOLD response corresponding to neural response of 1. The default
        is 0.05.
    fwhm_bold : float, optional
        BOLD point-spread width in mm. The default is 1.02.
    fwhm_noise : float, optional
        Measurement noise of BOLD response. The default is 0.001.
    a_mask : float, optional
        Major axis of elliptical mask. The default is 1000.
    b_mask : float, optional
        Minor axis of elliptical mask. The default is 1000.
    alpha_mask : float, optional
        Rotational angle of elliptical mask. The default is 0.
    path_white : str, optional
        Path to existing white noise image. The default is None.

    Returns
    -------
    arr_white : ndarray
        Initial white noise pattern.
    odc_img : ndarray
        Neural map.
    y_img : ndarray
        BOLD response with measurement noise.
    ymri_img : ndarray
        Sampled MRI signal.
    f_odc_fft : ndarray
        Anisotropic band-pass filter in spatial frequency representation.
    f_bold_fft : ndarray
        BOLD filter in spatial frequency representation.

    References
    -------
    .. [1] Rojer, AS, et al. Cat and monkey cortical columnar patterns modeled
    by bandpass-filtered 2D white noise, Biol Cybern 62(5), 381--391 (1990).
    .. [2] Chaimow, D, et al. Modeling and analysis of mechanisms underlying
    fMRI-based decoding of information conveyed in cortical columns, Neuroimage
    56(2), 627--642 (2011).

    """
    # get white noise
    if path_white:
        img = nb.load(path_white)
        arr_white = img.get_fdata()
    else:
        arr_white = get_white_2d(nx_sim, ny_sim, 0, 1)

    # get band-pass filter
    f_odc_fft = filter_odc_2d(nx_sim, ny_sim, fov_x, fov_y, rho, delta, epsilon, theta)

    # generate ODC pattern (neural map)
    odc_img = np.real(ifft2(fft2(arr_white) * f_odc_fft))
    odc_img = filter_sigmoid(odc_img, alpha)

    # BOLD response
    f_bold_fft = filter_bold_2d(nx_sim, ny_sim, fov_x, fov_y, fwhm_bold, beta)
    y_img = np.real(ifft2(fft2(odc_img) * f_bold_fft))

    # add measurement noise
    noise_img = get_white_2d(
        nx_sim, ny_sim, 0, fwhm_noise / (2 * np.sqrt(2 * np.log(2)))
    )
    y_img = noise_img + y_img

    # voxel sampling
    y_fft = fft2(y_img)

    kx_sample = np.round(nx_mri / 2).astype(int)
    ky_sample = np.round(ny_mri / 2).astype(int)

    ymri_fft = np.zeros((nx_mri, ny_mri), dtype=complex)
    ymri_fft[:kx_sample, :ky_sample] = y_fft[:kx_sample, :ky_sample]
    ymri_fft[:kx_sample, -1 : -ky_sample - 1 : -1] = y_fft[
        :kx_sample, -1 : -ky_sample - 1 : -1
    ]
    ymri_fft[-1 : -kx_sample - 1 : -1, :ky_sample] = y_fft[
        -1 : -kx_sample - 1 : -1, :ky_sample
    ]
    ymri_fft[-1 : -kx_sample - 1 : -1, -1 : -ky_sample - 1 : -1] = y_fft[
        -1 : -kx_sample - 1 : -1, -1 : -ky_sample - 1 : -1
    ]

    ymri_img = np.real(ifft2(ymri_fft))

    # mask voxel sampling
    ymri_img = ymri_img * mask_pattern_2d(
        np.shape(ymri_img)[0], np.shape(ymri_img)[1], a_mask, b_mask, alpha_mask
    )

    return arr_white, odc_img, y_img, ymri_img, f_odc_fft, f_bold_fft


def odc_1d(
    n_sim=1024,
    fov=20,
    n_mri=100,
    rho=0.5,
    delta=0.3,
    alpha=4,
    beta=0.05,
    fwhm_bold=1.02,
    fwhm_noise=0.001,
    a_mask=1000,
    b_mask=1000,
    path_white=None,
):
    """This function generates realistic ocular dominance patterns according to the
    model proposed by [1] in 1D. The rest follows similar to the ODC generation in 2D.

    Parameters
    ----------
    n_sim : int, optional
        Array size of the simulated patch in x-direction. The default is 1024.
    fov : float, optional
        Field of view in x-direction (mm). The default is 20.
    n_mri : int, optional
        Array size of the MR image in x-direction. The default is 100.
    rho : float, optional
        Main spatial frequency determining columnar width in cycles/mm. The
        default is 0.5.
    delta : float, optional
        Variations orthogonal to ODC bands in cycles/mm (irregularity). The
        default is 0.3.
    alpha : float, optional
        Sharpness parameter of the sigmoidal filter. The default is 4.
    beta : float, optional
        Maximal BOLD response corresponding to neural response of 1. The default
        is 0.05.
    fwhm_bold : float, optional
        BOLD point-spread width in mm. The default is 1.02.
    fwhm_noise : float, optional
        Measurement noise of BOLD response. The default is 0.001.
    a_mask : float, optional
        Left side length of mask. The default is 1000.
    b_mask : float, optional
        Right side length of mask. The default is 1000.
    path_white : str, optional
        Path to existing white noise image. The default is None.

    Returns
    -------
    arr_white : ndarray
        Initial white noise pattern.
    odc_img : ndarray
        Neural map.
    y_img : ndarray
        BOLD response with measurement noise.
    ymri_img : ndarray
        Sampled MRI signal.
    f_odc_fft : ndarray
        Anisotropic band-pass filter in spatial frequency representation.
    f_bold_fft : ndarray
        BOLD filter in spatial frequency representation.

    References
    -------
    .. [1] Rojer, AS, et al. Cat and monkey cortical columnar patterns modeled
    by bandpass-filtered 2D white noise, Biol Cybern 62(5), 381--391 (1990).

    """
    # add path of the executed script to the interpreter's search path
    sys.path.append(sys.argv[0])

    # get white noise
    if path_white:
        img_white = nb.load(path_white)
        arr_white = img_white.get_fdata()
    else:
        arr_white = get_white_1d(n_sim, 0, 1)

    # get band-pass filter
    f_odc_fft = filter_odc_1d(n_sim, fov, rho, delta)

    # generate ODC pattern (neural map)
    odc_img = np.real(ifft(fft(arr_white) * f_odc_fft))
    odc_img = filter_sigmoid(odc_img, alpha)

    # BOLD response
    f_bold_fft = filter_bold_1d(n_sim, fov, fwhm_bold, beta)
    y_img = np.real(ifft(fft(odc_img) * f_bold_fft))

    # add measurement noise
    noise_img = get_white_1d(n_sim, 0, fwhm_noise / (2 * np.sqrt(2 * np.log(2))))
    y_img = noise_img + y_img

    # voxel sampling
    y_fft = fft(y_img)

    k_sample = np.round(n_mri / 2).astype(int)

    ymri_fft = np.zeros(n_mri, dtype=complex)
    ymri_fft[:k_sample] = y_fft[:k_sample]
    ymri_fft[-1 : -k_sample - 1 : -1] = y_fft[-1 : -k_sample - 1 : -1]

    ymri_img = np.real(ifft(ymri_fft))

    # mask voxel sampling
    ymri_img = ymri_img * mask_pattern_1d(len(ymri_img), a_mask, b_mask)

    return arr_white, odc_img, y_img, ymri_img, f_odc_fft, f_bold_fft


def pattern_2d(
    nx_sim=1024,
    ny_sim=1024,
    fov_x=20,
    fov_y=20,
    nx_mri=100,
    ny_mri=100,
    omega_x=0.5,
    omega_y=0.5,
    phi_x=0,
    phi_y=0,
    theta=0,
    rect_shape=True,
    beta=0.05,
    fwhm_bold=1.02,
    fwhm_noise=0.001,
    a_mask=1000,
    b_mask=1000,
    alpha_mask=0,
):
    """This function geometrical stripe pattern in 2D. The rest is similar to the ODC
    pattern generation.

    Parameters
    ----------
    nx_sim : int, optional
        Array size of the simulated patch in x-direction. The default is 1024.
    ny_sim : int, optional
        Array size of the simulated patch in y-direction. The default is 1024.
    fov_x : float, optional
        Field of view in x-direction (mm). The default is 20.
    fov_y : float, optional
        Field of view in y-direction (mm). The default is 20.
    nx_mri : int, optional
        Array size of the MR image in x-direction. The default is 100.
    ny_mri : int, optional
        Array size of the MR image in y-direction. The default is 100.
    omega_x : float, optional
        Frequency in x-direction in cycles/mm. The default is 0.5.
    omega_y : float, optional
        Frequency in y-direction in cycles/mm. The default is 0.5.
    phi_x : float, optional
        Phase in x-direction in deg. The default is 0.
    phi_y : float, optional
        Phase in y-direction in deg. The default is 0.
    theta : float, optional
        Rotation angle in deg. The default is 0.
    rect_shape : bool, optional
        Rectangular or sine pattern. The default is True.
    beta : float, optional
        Maximal BOLD response corresponding to neural response of 1. The default
        is 0.05.
    fwhm_bold : float, optional
        BOLD point-spread width in mm. The default is 1.02.
    fwhm_noise : float, optional
        Measurement noise of BOLD response. The default is 0.001.
    a_mask : float, optional
        Major axis of elliptical mask. The default is 1000.
    b_mask : float, optional
        Minor axis of elliptical mask. The default is 1000.
    alpha_mask : float, optional
        Rotational angle of elliptical mask. The default is 0.

    Returns
    -------
    arr_pattern : ndarray
        Neural map.
    y_img : ndarray
        BOLD response with measurement noise.
    ymri_img : ndarray
        Sampled MRI signal.
    f_bold_fft : ndarray
        BOLD filter in spatial frequency representation.

    """
    # generate geometrical pattern
    x = np.linspace(0, fov_x, nx_sim)
    y = np.linspace(0, fov_y, ny_sim)

    x_grid, y_grid = np.meshgrid(x, y)

    # compute rotation
    theta = np.radians(theta)
    c, s = np.cos(theta), np.sin(theta)
    x_rot = c * x_grid - s * y_grid
    y_rot = s * x_grid + c * y_grid

    # convert phase from deg to rad
    phi_x = np.radians(phi_x)
    phi_y = np.radians(phi_y)

    cos_x = np.cos(2 * np.pi * omega_x * x_rot + phi_x)
    cos_y = np.cos(2 * np.pi * omega_y * y_rot + phi_y)

    # merge pattern of x- and y-direction
    arr_pattern = cos_x * cos_y

    # generate rect or sine pattern
    if rect_shape:
        arr_pattern[arr_pattern > 0] = 1
        arr_pattern[arr_pattern != 1] = -1

    # BOLD response
    f_bold_fft = filter_bold_2d(nx_sim, ny_sim, fov_x, fov_y, fwhm_bold, beta)
    y_img = np.real(ifft2(fft2(arr_pattern) * f_bold_fft))

    # add measurement noise
    noise_img = get_white_2d(
        nx_sim, ny_sim, 0, fwhm_noise / (2 * np.sqrt(2 * np.log(2)))
    )
    y_img = noise_img + y_img

    # voxel sampling
    y_fft = fft2(y_img)

    kx_sample = np.round(nx_mri / 2).astype(int)
    ky_sample = np.round(ny_mri / 2).astype(int)

    ymri_fft = np.zeros((nx_mri, ny_mri), dtype=complex)
    ymri_fft[:kx_sample, :ky_sample] = y_fft[:kx_sample, :ky_sample]
    ymri_fft[:kx_sample, -1 : -ky_sample - 1 : -1] = y_fft[
        :kx_sample, -1 : -ky_sample - 1 : -1
    ]
    ymri_fft[-1 : -kx_sample - 1 : -1, :ky_sample] = y_fft[
        -1 : -kx_sample - 1 : -1, :ky_sample
    ]
    ymri_fft[-1 : -kx_sample - 1 : -1, -1 : -ky_sample - 1 : -1] = y_fft[
        -1 : -kx_sample - 1 : -1, -1 : -ky_sample - 1 : -1
    ]

    ymri_img = np.real(ifft2(ymri_fft))

    # mask voxel sampling
    ymri_img = ymri_img * mask_pattern_2d(
        np.shape(ymri_img)[0], np.shape(ymri_img)[1], a_mask, b_mask, alpha_mask
    )

    return arr_pattern, y_img, ymri_img, f_bold_fft


def pattern_1d(
    n_sim=1024,
    fov=20,
    n_mri=100,
    omega=0.5,
    phi=0,
    rect_shape=True,
    beta=0.05,
    fwhm_bold=1.02,
    fwhm_noise=0.001,
    a_mask=1000,
    b_mask=1000,
):
    """This function geometrical stripe pattern in 1D. The rest is similar to the stripe
    pattern generation in 2D.

    Parameters
    ----------
    n_sim : int, optional
        Array size of the simulated patch. The default is 1024.
    fov : float, optional
        Field of view (mm). The default is 20.
    n_mri : int, optional
        Array size of the MR image. The default is 100.
    omega : float, optional
        Frequency in cycles/mm. The default is 0.5.
    phi : float, optional
        Phase in deg. The default is 0.
    rect_shape : bool, optional
        Rectangular or sine pattern. The default is True.
    beta : float, optional
        Maximal BOLD response corresponding to neural response of 1. The default
        is 0.05.
    fwhm_bold : float, optional
        BOLD point-spread width in mm. The default is 1.02.
    fwhm_noise : float, optional
        Measurement noise of BOLD response. The default is 0.001.
    a_mask : float, optional
        Left side length of mask. The default is 1000.
    b_mask : float, optional
        Right side length of mask. The default is 1000.

    Returns
    -------
    arr_pattern : ndarray
        Neural map.
    y_img : ndarray
        BOLD response with measurement noise.
    ymri_img : ndarray
        Sampled MRI signal.
    f_bold_fft : ndarray
        BOLD filter in spatial frequency representation.

    """
    # add path of the executed script to the interpreter's search path
    sys.path.append(sys.argv[0])

    # generate geometrical pattern
    x = np.linspace(0, fov, n_sim)

    # convert phase from deg to rad
    phi = np.radians(phi)

    arr_pattern = np.cos(2 * np.pi * omega * x + phi)

    # generate rect or sine pattern
    if rect_shape:
        arr_pattern[arr_pattern > 0] = 1
        arr_pattern[arr_pattern != 1] = -1

    # BOLD response
    f_bold_fft = filter_bold_1d(n_sim, fov, fwhm_bold, beta)
    y_img = np.real(ifft(fft(arr_pattern) * f_bold_fft))

    # add measurement noise
    noise_img = get_white_1d(n_sim, 0, fwhm_noise / (2 * np.sqrt(2 * np.log(2))))
    y_img = noise_img + y_img

    # voxel sampling
    y_fft = fft(y_img)

    k_sample = np.round(n_mri / 2).astype(int)

    ymri_fft = np.zeros(n_mri, dtype=complex)

    ymri_fft[:k_sample] = y_fft[:k_sample]
    ymri_fft[-1 : -k_sample - 1 : -1] = y_fft[-1 : -k_sample - 1 : -1]

    ymri_img = np.real(ifft(ymri_fft))

    # mask voxel sampling
    ymri_img = ymri_img * mask_pattern_1d(len(ymri_img), a_mask, b_mask)

    return arr_pattern, y_img, ymri_img, f_bold_fft


def filter_bold_2d(nx, ny, fovx, fovy, fwhm, beta):
    """This defines the BOLD blurring in 2D. It is defined as modulation transfer
    function in spatial frequency space. The filter is the Fourier transform of a
    two-dimensional Gaussian kernel in image space.

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
    fwhm : float
        Full width at half maximum of the Gaussian kernel.
    beta : float
        Maximum BOLD response for neural response of 1.

    Returns
    -------
    f : ndarray
        Band-pass filter array in spatial frequency space.

    """
    # get maximum k-space coordinate in x- and y-direction
    kx_max = nx / (2 * fovx)
    ky_max = ny / (2 * fovy)

    # get k-space axes
    kx = np.linspace(-kx_max, kx_max, nx)
    ky = np.linspace(-ky_max, ky_max, ny)

    # define two-dimensional k-space grid
    kx_grid, ky_grid = np.meshgrid(kx, ky)

    # convert to polar coordinates
    k_r = np.sqrt(kx_grid**2 + ky_grid**2)

    # BOLD modulation transfer function
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))

    f = beta * np.exp(-2 * np.pi**2 * sigma**2 * k_r**2)

    # shift zero k-space line to the left border to match the numpy fft
    # convention
    f = fftshift(f)

    return f


def filter_bold_1d(n, fov, fwhm, beta):
    """This defines the BOLD blurring in 1D. It is defined as modulation transfer
    function in spatial frequency space. The filter is the Fourier transform of a
    one-dimensional Gaussian kernel in image space.

    Parameters
    ----------
    n : int
        Array size.
    fov : float
        Field of view (mm).
    fwhm : float
        Full width at half maximum of the Gaussian kernel.
    beta : float
        Maximum BOLD response for neural response of 1.

    Returns
    -------
    f : ndarray
        Band-pass filter array in spatial frequency space.

    """
    # get maximum k-space coordinate in x- and y-direction
    k_max = n / (2 * fov)

    # get k-space axes
    k = np.linspace(-k_max, k_max, n)

    # BOLD modulation transfer function
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))

    f = beta * np.exp(-2 * np.pi**2 * sigma**2 * k**2)

    # shift zero k-space line to the left border to match the numpy fft
    # convention
    f = fftshift(f)

    return f


def filter_odc_2d(nx, ny, fovx, fovy, rho, delta, epsilon, theta):
    """This function defines the 2D ocular dominance column filter (anisotropic
    band-pass filter) in spatial frequency space which is taken from [1]. The filter is
    normalised in order to have the same variance as the white noise source image.

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
    """This function defines the 1D ocular dominance column filter (anisotropic
    band-pass filter) in spatial frequency space which is taken from [1]. The filter is
    normalised in order to have the same variance as the white noise source image.

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


def filter_sigmoid(x, alpha=4):
    """Nonlinear sigmoid filter to adjust sharpness of ocular dominance. If sharpness
    factor is set to zero, nothing is changed on the input image. The sigmoidal filter
    outputs an array with intensity range [-1,1].

    Parameters
    ----------
    x : ndarray
        Input image.
    alpha : float, optional
        Sharpness factor. The default is 4.

    Returns
    -------
    ndarray
        Filtered image.

    """
    if alpha:
        # return 1 / (1 + np.exp(-alpha*x))
        return 2 / (1 + np.exp(-alpha * x)) - 1
    else:
        return x


def get_white_2d(nx, ny, mu, sigma):
    """Creates a 2D array with Gaussian noise.

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
    """Creates a 1D array with Gaussian noise.

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


def mask_pattern_2d(nx, ny, a, b, alpha):
    """This function computes a mask to occlude parts of the columnar pattern. The mask
    has elliptical shape and can be adjusted by its major and minor axes and a
    rotational angle.

    Parameters
    ----------
    nx : int
        Array size in x-direction.
    ny : int
        Array size in y-direction.
    a : float
        Major axis of ellipse.
    b : float
        Minor axis of ellipse.
    alpha : float
        Rotation in degrees.

    Returns
    -------
    mask : ndarray
        Binary array containing mask.

    """
    # compute meshgrid
    x = np.linspace(-nx / 2, nx / 2, nx)
    y = np.linspace(-ny / 2, ny / 2, ny)
    x_mesh, y_mesh = np.meshgrid(x, y)

    # rotate meshgrid
    x_rot = x_mesh * np.cos(math.radians(alpha)) - y_mesh * np.sin(math.radians(alpha))
    y_rot = x_mesh * np.sin(math.radians(alpha)) + y_mesh * np.cos(math.radians(alpha))

    # compute ellipse
    mask = x_rot**2 / a**2 + y_rot**2 / b**2

    # mask ellipse
    mask[mask > 1] = 0
    mask[mask != 0] = 1

    return mask


def mask_pattern_1d(n, a, b):
    """This function computes a mask to occlude parts of the columnar pattern.

    Parameters
    ----------
    n : int
        Array size.
    a : float
        Mask left side.
    b : float
        Mask right side.

    Returns
    -------
    mask : ndarray
        Binary array containing mask.

    """
    # compute meshgrid
    mask = np.linspace(-n / 2, n / 2, n)

    # mask ellipse
    mask[mask < -a] = 0
    mask[mask > b] = 0
    mask[mask != 0] = 1

    return mask


def regrid_zero_2d(data_array, nx_new, ny_new):
    """This function loads a two-dimensional numpy array and performs spatial
    interpolation using zero padding in spatial frequency space. If the input array has
    an odd array size, the padding size on the right size is decreased by one to match
    the output array size and leave the k-space center line at the same place as for
    even array sizes.

    Parameters
    ----------
    data_array : ndarray
        Input array.
    nx_new : int
        Resolution of output array in x-direction (even number).
    ny_new : int
        Resolution of output array in y-direction (even number).

    Returns
    -------
    data_array_new : ndarray
        Spatial interpolation of input array.

    """
    # size of input array
    data_size = np.shape(data_array)

    nx_old = data_size[0]
    ny_old = data_size[1]

    # zero padding size
    padx_left = int(nx_new / 2) - int(nx_old / 2)
    pady_left = int(ny_new / 2) - int(ny_old / 2)
    padx_right = int(nx_new / 2) - int(nx_old / 2)
    pady_right = int(ny_new / 2) - int(ny_old / 2)

    # check if old data arrays is even/odd
    if nx_old % 2:
        padx_right -= 1

    if ny_old % 2:
        pady_right -= 1

    # Fourier transform of input data
    data_array_fft = fft2(data_array)
    data_array_fft = fftshift(data_array_fft)

    # pad zeros in spatial frequency space
    data_array_fft = np.pad(
        data_array_fft, ((padx_left, padx_right), (pady_left, pady_right)), "constant"
    )

    data_array_fft = fftshift(data_array_fft)

    # inverse Fourier transform of zero padded spatial frequency representation
    data_array_new = ifft2(data_array_fft)
    data_array_new = np.real(data_array_new)

    return data_array_new


def regrid_zero_1d(data_array, n_new):
    """This function loads a one-dimensional numpy array and performs spatial
    interpolation using zero padding spatial frequency space. If the input array has an
    odd array size, the padding size on the right size is decreased by one to match the
    output array size and leave the k-space center line at the same place as for even
    array sizes.

    Parameters
    ----------
    data_array : ndarray
        Input array.
    n_new : int
        Resolution of output array (even number).

    Returns
    -------
    data_array_new : ndarray
        Spatial interpolation of input array.

    """
    # size of input array
    data_size = np.shape(data_array)

    n_old = data_size[0]

    # zero padding size
    pad_left = int(n_new / 2) - int(n_old / 2)
    pad_right = int(n_new / 2) - int(n_old / 2)

    # check if old or new data arrays are even/odd
    if n_old % 2:
        pad_right -= 1

    # Fourier transform of input data
    data_array_fft = fft(data_array)
    data_array_fft = fftshift(data_array_fft)

    # pad zeros in spatial frequency space
    data_array_fft = np.pad(data_array_fft, (pad_left, pad_right), "constant")

    data_array_fft = fftshift(data_array_fft)

    # inverse Fourier transform of zero padded spatial frequency representation
    data_array_new = ifft(data_array_fft)
    data_array_new = np.real(data_array_new)

    return data_array_new


def regrid_2d(data_array, nx, ny):
    """This function loads a two-dimensional numpy array and interpolates it to a new
    grid array with defined sizes using nearest neighbour interpolation.

    Parameters
    ----------
    data_array : ndarray
        Input array.
    nx : int
        Resolution of output array in x-direction.
    ny : int
        Resolution of output array in y-direction.

    Returns
    -------
    data_array_new : ndarray
        Output array.

    """
    # size of input array
    data_size = np.shape(data_array)

    # grid of input array
    x = np.arange(0, data_size[0])
    y = np.arange(0, data_size[1])
    xgrid, ygrid = np.meshgrid(x, y)
    xgrid = np.reshape(xgrid, np.size(xgrid))
    ygrid = np.reshape(ygrid, np.size(ygrid))
    xi_old = np.stack((xgrid, ygrid), 1)

    # grid of output array
    x_new = np.linspace(0, data_size[0], nx)
    y_new = np.linspace(0, data_size[1], ny)
    xgrid_new, ygrid_new = np.meshgrid(x_new, y_new)
    xgrid_new = np.reshape(xgrid_new, np.size(xgrid_new))
    ygrid_new = np.reshape(ygrid_new, np.size(ygrid_new))
    # xi_new = np.floor(np.stack((xgrid_new,ygrid_new),1))
    xi_new = np.stack((xgrid_new, ygrid_new), 1)

    # values of input array
    data_array = np.reshape(data_array, np.size(data_array))

    # grid values from old grid to new grid
    data_array_new = griddata(xi_old, data_array, xi_new, method="nearest")
    data_array_new = np.reshape(data_array_new, (nx, ny))

    return data_array_new


def regrid_1d(data_array, n):
    """This function loads a one-dimensional numpy array and interpolates it to a new
    grid array with defined sizes using nearest neighbour interpolation.

    Parameters
    ----------
    data_array : ndarray
        Input array.
    n : int
        Resolution of output array.

    Returns
    -------
    data_array_new : ndarray
        Array values interpolated to new grid.

    """
    # size of input array
    data_size = np.shape(data_array)

    # grid of input array
    xi_old = np.arange(0, data_size[0])

    # grid of output array
    # xi_new = np.floor(np.linspace(0,data_size[0],N))
    xi_new = np.linspace(0, data_size[0], n)

    # grid values from old grid to new grid
    data_array_new = griddata(xi_old, data_array, xi_new, method="nearest")

    return data_array_new


def pattern_corr(data_array1, data_array2):
    """This function calculates the pearson correlation between two arrays of the same
    size.

    Parameters
    ----------
    data_array1 : ndarray
        Input array1.
    data_array2 : ndarray
        Input array2.

    Returns
    -------
    tuple
        Pearson coefficient

    """
    # reshape data arrays into a vector
    x = np.reshape(data_array1, (np.size(data_array1), 1))
    y = np.reshape(data_array2, (np.size(data_array2), 1))

    return pearsonr(x, y)
