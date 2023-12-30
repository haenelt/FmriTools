# -*- coding: utf-8 -*-

import sys

import nibabel as nb
import numpy as np
from numpy.fft import fft, fft2, ifft, ifft2

from ..simulation.filter_bold import filter_bold_1d, filter_bold_2d
from ..simulation.filter_odc import filter_odc_1d, filter_odc_2d
from ..simulation.filter_sigmoid import filter_sigmoid
from ..simulation.get_white import get_white_1d, get_white_2d
from ..simulation.mask_pattern import mask_pattern_1d, mask_pattern_2d


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
    """ODC 2D.

    This function generates realistic ocular dominance patterns according to the
    model proposed by [1] in 2D. Implementation follows [2] and most of the
    default values for the columnar pattern are taken from this publication.
    First, a white noise pattern is defined and an anisotropic band-pass filter
    and a non-linear sigmoidal filter are applied. This neural map is converted
    to a BOLD response map by applying a Gaussian filter. White noise is added
    to simulate uncorrelated measurement noise to the signal. The MRI sampling
    procedure is considered by only inner k-space lines of the spatial frequency
    representation of the BOLD reseponse map.

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
    """ODC 1D.

    This function generates realistic ocular dominance patterns according to the
    model proposed by [1] in 1D. The rest follows similar to the ODC generation
    in 2D.

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
