# -*- coding: utf-8 -*-

import numpy as np
from numpy.fft import fftshift


def get_bandpass_filter(
    nx,
    ny,
    fovx,
    fovy,
    kcut_low=0,
    kcut_high=1,
    apply_fftshift=False,
    k_fwhm=0,
    theta_fwhm=0,
    theta1=0,
    theta2=180,
):
    """Get bandpass filter.

    This function generates a bandpass filter for an input image defined by
    lower and upper cutoff frequencies. Additionally, Gaussian attenuation of
    border can be appled and the filter can be only defined to a specific
    k-space region defined by lower and upper angle thresholds.

    Parameters
    ----------
    nx : int
        Matrix size of input array in x-direction.
    ny : int
        Matrix size of input array in y-direction.
    fovx : float
        Field of view of input array in x-direction in mm.
    fovy : float
        Field of view of input array in y-direction in mm.
    kcut_low : float, optional
        Lower spatial cutoff frequency in cycles/mm. The default is 0.
    kcut_high : föpat, optional
        Upper spatial cutoff frequency in cycles/mm. The default is 1.
    apply_fftshift : bool, optional
        FFTshift to stay in the convention of spatial frequencies in numpy
        arrays. The default is False.
    k_fwhm : float, optional
        Full-width at half maximum of gaussian filter in frequency direction.
        The default is 0.
    theta_fwhm : float, optional
        Full-width at half maximum of gaussian filter in angle direction. The
        default is 0.
    theta1 : float, optional
        Lower cutoff angle in deg [0,180]. The default is 0.
    theta2 : float, optional
        Higher cutoff angle in deg [0,180]. The default is 180.

    Returns
    -------
    B : ndarray
        Spatial frequency filter.

    """

    # parameters of gaussian
    beta = 1
    k_sigma = k_fwhm / (2 * np.sqrt(2 * np.log(2)))
    theta_sigma = theta_fwhm / (2 * np.sqrt(2 * np.log(2)))

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

    k_pol = np.arctan2(ky_grid, kx_grid)
    k_pol[k_pol < 0] = k_pol[k_pol < 0] + np.pi
    k_pol = k_pol / np.pi * 180

    # frequency filter
    b_r = k_r.copy()
    b_r[np.where(np.logical_and(b_r >= kcut_low, b_r <= kcut_high))] = np.nan
    b_r[~np.isnan(b_r)] = 0
    b_r[b_r != 0] = 1

    # angle filter
    b_pol = k_pol.copy()
    if theta2 > theta1:
        b_pol[np.where(np.logical_and(b_pol >= theta1, b_pol <= theta2))] = np.nan
        b_pol[~np.isnan(b_pol)] = 0
        b_pol[b_pol != 0] = 1
    else:
        b_pol[np.where(np.logical_and(b_pol <= theta1, b_pol >= theta2))] = np.nan
        b_pol[np.isnan(b_pol)] = 0
        b_pol[b_pol != 0] = 1
        b_pol[k_pol == 0] = 1  # important to also fill the phase gap

    # filter edges
    if k_fwhm != 0:
        b1 = (
            beta
            / (np.sqrt(2 * np.pi) * k_sigma)
            * np.exp(-((k_r - kcut_low) ** 2) / (2 * k_sigma**2))
        )
        b2 = (
            beta
            / (np.sqrt(2 * np.pi) * k_sigma)
            * np.exp(-((k_r - kcut_high) ** 2) / (2 * k_sigma**2))
        )

        b_temp = b1 + b2
        b_temp2 = b_temp.copy()
        b_temp2[b_r == 1] = 0
        b_max = np.max(b_temp2)
        b_temp[b_r == 1] = b_max
        b_r = b_temp.copy()

        # normalize filter
        b_r = (b_r - np.min(b_r)) / (np.max(b_r) - np.min(b_r))

    if theta_fwhm != 0 and theta2 - theta1 < 180:
        angle1 = np.mod(k_pol - theta1, 180)
        angle2 = np.mod(theta1 - k_pol, 180)
        angle3 = np.mod(k_pol - theta2, 180)
        angle4 = np.mod(theta2 - k_pol, 180)

        b_pol1a = (
            beta
            / (np.sqrt(2 * np.pi) * theta_sigma)
            * np.exp(-(angle1**2) / (2 * theta_sigma**2))
        )
        b_pol1b = (
            beta
            / (np.sqrt(2 * np.pi) * theta_sigma)
            * np.exp(-(angle2**2) / (2 * theta_sigma**2))
        )
        b_pol2a = (
            beta
            / (np.sqrt(2 * np.pi) * theta_sigma)
            * np.exp(-(angle3**2) / (2 * theta_sigma**2))
        )
        b_pol2b = (
            beta
            / (np.sqrt(2 * np.pi) * theta_sigma)
            * np.exp(-(angle4**2) / (2 * theta_sigma**2))
        )

        b_temp = b_pol1a + b_pol1b + b_pol2a + b_pol2b
        b_temp2 = b_temp.copy()
        b_temp2[b_pol == 1] = 0
        b_max = np.max(b_temp2)
        b_temp[b_pol == 1] = b_max
        b_pol = b_temp.copy()

        # normalize filter
        b_pol = (b_pol - np.min(b_pol)) / (np.max(b_pol) - np.min(b_pol))

    # compose frequency and angle filter
    b = b_r * b_pol

    # shift zero k-space line to the left border to match the numpy fft
    # convention
    if apply_fftshift is True:
        b = fftshift(b)

    return b
