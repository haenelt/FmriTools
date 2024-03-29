# -*- coding: utf-8 -*-
"""Spatial properties of image arrays."""

import copy
import os
import random
import re

import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
import pydicom
from nibabel.freesurfer.io import read_label
from numpy.fft import fft, fft2, fftshift
from numpy.random import shuffle
from scipy.signal import find_peaks
from scipy.stats import levene, ttest_ind

from ..io.filename import get_filename
from ..utils.metrics import calc_acorr

__all__ = [
    "analyze_acorr",
    "analyze_alff_between_conditions",
    "analyze_alff_between_stripes",
    "analyze_fft",
    "get_pca",
    "get_fft",
    "get_bandpass",
    "get_rf_pulse_bw",
    "get_bandwidth_per_px",
]


def analyze_acorr(arr, fovx, fovy, xv, yv, p_min=0.01, p_max=0.5, nsample=1000):
    """This function computes the normalized autocorrelation (NAC) from a 2D input array
    and estimates the width of the central peak and the distance to its first neighbor
    peak along a defined projection lines sampled with nearest neighbor interpolation.
    The width of the central peak is defined as the width at zero point.

    Parameters
    ----------
    arr : ndarray
        2D input array.
    fovx : float
        Field of view in x-direction (mm).
    fovy : float
        Field of view in y-firection (mm).
    xv : float
        x-coordinate of pca eigenvector.
    yv : float
        y-coordinate of pca eigenvector.
    p_min : float, optional
        Minimum prominence for peak detection. The default is 0.01.
    p_max : float, optional
        Maximum prominence for peak detection. The default is 0.5.
    nsample : int, optional
        Number of sampling points along projection line. The default is 1000.

    Returns
    -------
    fwhm_central : float
        FWHM of NAC central peak in mm along projection axis.
    d_neighbor : float
        Distance to first neighbor peak in mm along projection axis.
    p_neighbor : float
        Power of first neighbor peak along projection axis.
    d : float
        Lag in mm along projection axis.
    acorr_line : ndarray
        NAC along projection axis.

    """
    # add one if nsample is an odd integer
    if np.mod(nsample, 2) > 0:
        nsample += 1

    # get size of array
    x_size = np.shape(arr)[0]
    y_size = np.shape(arr)[1]

    # get coordinates of array
    x = (fovx / 2) * np.linspace(-1, 1, x_size)
    y = (fovy / 2) * np.linspace(-1, 1, y_size)
    y_mesh, x_mesh = np.meshgrid(y, x)

    # get projection line coordinates from pca eigenvector. Because we
    # interpolate using nearest neighbors the axis enc point is <size>-1.
    if np.abs(xv) < np.abs(yv):
        y_line = np.linspace(0, y_size - 1, nsample)
        x_line = xv / yv * y_line + y_size / 2
        x_line = x_line + x_size / 2 - x_line[int(nsample / 2)]
    else:
        x_line = np.linspace(0, x_size - 1, nsample)
        y_line = yv / xv * x_line + x_size / 2
        y_line = y_line + y_size / 2 - y_line[int(nsample / 2)]

    # check that no line steps over array border
    if np.any(y_line > y_size - 1):
        x_line[y_line > y_size - 1] = np.nan
        y_line[y_line > y_size - 1] = np.nan
        x_line = x_line[~np.isnan(x_line)]
        y_line = y_line[~np.isnan(y_line)]

    if np.any(y_line < 0):
        x_line[y_line < 0] = np.nan
        y_line[y_line < 0] = np.nan
        x_line = x_line[~np.isnan(x_line)]
        y_line = y_line[~np.isnan(y_line)]

    if np.any(x_line > x_size - 1):
        y_line[x_line > x_size - 1] = np.nan
        x_line[x_line > x_size - 1] = np.nan
        x_line = x_line[~np.isnan(x_line)]
        y_line = y_line[~np.isnan(y_line)]

    if np.any(y_line < 0):
        y_line[x_line < 0] = np.nan
        x_line[x_line < 0] = np.nan
        x_line = x_line[~np.isnan(x_line)]
        y_line = y_line[~np.isnan(y_line)]

    xx_line = x_mesh[np.round(x_line).astype(int), np.round(y_line).astype(int)]
    yy_line = y_mesh[np.round(x_line).astype(int), np.round(y_line).astype(int)]

    # get final distances
    d = np.sqrt(xx_line**2 + yy_line**2)

    # get minimum to shift mid point to origin
    d_min = np.argwhere(np.min(d) == d)

    if np.size(d_min) > 1:
        d_min = int(d_min[-1])
    else:
        d_min = int(d_min)

    # all coordinates left from origin get negative signs
    d[:d_min] = -d[:d_min]

    # get autocorrelation
    array_acorr = calc_acorr(arr)
    acorr_line = array_acorr[np.round(x_line).astype(int), np.round(y_line).astype(int)]

    # fwhm
    acorr_line_max = np.argwhere(np.max(acorr_line) == acorr_line)
    if np.size(acorr_line_max) > 1:
        acorr_line_max = int(acorr_line_max[0])
    else:
        acorr_line_max = int(acorr_line_max)

    # FWHM will be defined at half maximum. N.B., this underestimates the
    # columnar width in case of pure sinusoidal oscillation where the width
    # would be determined by taking the FWHM at zero
    acorr_line_middle = copy.deepcopy(acorr_line_max)
    while True:
        acorr_line_middle += 1
        if acorr_line_middle > len(acorr_line) - 1:
            acorr_line_middle = np.nan
            break
        elif acorr_line[acorr_line_middle] < 0.5:
            break

    # compute fwhm if middle point is found
    if ~np.isnan(acorr_line_middle):
        fwhm_central = 2 * np.abs(d[acorr_line_max] - d[acorr_line_middle])
    else:
        fwhm_central = np.nan

    # spacing to neighbor
    peak = find_peaks(acorr_line, prominence=(p_min, p_max))[0]

    if len(peak) <= 1:
        d_neighbor = np.nan
        p_neighbor = np.nan
    else:
        peak_temp = peak - acorr_line_max
        peak_temp = np.abs(peak_temp)
        peak_temp = peak_temp.astype(float)
        peak_temp[peak_temp == 0] = np.NaN
        peak_temp = int(np.argwhere(np.nanmin(peak_temp) == peak_temp)[0])

        p_neighbor = acorr_line[peak[peak_temp]]
        d_neighbor = np.abs(d[acorr_line_max] - d[peak[peak_temp]])

    return fwhm_central, d_neighbor, p_neighbor, d, acorr_line


def analyze_alff_between_conditions(
    input_label, input_contrast1, input_contrast2, input_rest, min_contrast, nvert
):
    """Comparison of resting-state between different stripes populations. Mask stripes
    from an input contrast within a label ROI and randomly select <nvert> vertices
    within the final mask. An independent samples t-test is computed (or Welch's test if
    Levene's test is significant).

    Parameters
    ----------
    input_label : str
        Input label file to define the region of interest.
    input_contrast1 : str
        First input contrast data for masking.
    input_contrast2 : str
         Second input contrast data for masking.
    input_rest : str
        Input resting-state data (alff or falff).
    min_contrast : float
        Minimum resting-state fluctuation within the mask.
    nvert : int
        Number of selected vertices.

    Returns
    -------
    rest1 : ndarray
        Resting-state data within condition1.
    rest2 : ndarray
        Resting-state data within condition2.
    t : float
        t-score from independent samples t-test.
    p : float
        p-value from independent samples t-test.
    p_levene : float
        p-value from Levene's test.

    """
    # load data
    label, _ = read_label(input_label, read_scalars=True)
    contrast1 = np.squeeze(nb.load(input_contrast1).get_fdata())
    contrast2 = np.squeeze(nb.load(input_contrast2).get_fdata())
    rest = np.squeeze(nb.load(input_rest).get_fdata())

    # mask data
    contrast1[contrast1 < min_contrast] = np.NaN
    contrast2[contrast2 < min_contrast] = np.NaN

    contrast1[~label] = np.NaN
    contrast2[~label] = np.NaN
    contrast1[~np.isnan(contrast1)] = 1
    contrast2[~np.isnan(contrast2)] = 1

    # get resting-state data in mask
    rest1 = rest * contrast1
    rest1 = rest1[~np.isnan(rest1)]
    rest1 = rest1[rest1 != np.min(rest1)]

    rest2 = rest * contrast2
    rest2 = rest2[~np.isnan(rest2)]
    rest2 = rest2[rest2 != np.min(rest2)]

    # select random number of vertices
    label_shuffled1 = random.sample(range(0, len(rest1)), nvert)
    label_shuffled2 = random.sample(range(0, len(rest2)), nvert)

    rest1 = rest1[label_shuffled1]
    rest2 = rest2[label_shuffled2]

    # independent samples t-test
    # Levene's test is run to check for equal variances. If variances are not
    # equal, Welch's t-test is performed.
    _, p_levene = levene(rest1, rest2)
    if p_levene < 0.05:
        t, p = ttest_ind(rest1, rest2, equal_var=False)
    else:
        t, p = ttest_ind(rest1, rest2, equal_var=True)

    return rest1, rest2, t, p, p_levene


def analyze_alff_between_stripes(
    input_label, input_contrast, input_rest, min_contrast, nvert
):
    """Comparison of resting-state data within and between V2 stripes. Mask stripes
    within a ROI from a label with a thresholded contrast and select randomly <nvert>
    vertices either within or between stripes. An independent samples t-test is computed
    (or Welch's test if Levene's test is significant).

    Parameters
    ----------
    input_label : str
         Input label file to define the region of interest.
    input_contrast : str
        Input contrast data for masking.
    input_rest : str
        Input resting-state data (alff or falff).
    min_contrast : float
        Minimum contrast for masking (t-score).
    nvert : int
        Number of selected vertices.

    Returns
    -------
    rest_pos : ndarray
        Resting-state data within stripes.
    rest_neg : ndarray
        Resting-state data between stripes.
    t : float
        t-score from independent samples t-test.
    p : float
        p-value from independent samples t-test.
    p_levene : float
        p-value from Levene's test.

    """
    # load data
    label, _ = read_label(input_label, read_scalars=True)
    contrast = np.squeeze(nb.load(input_contrast).get_fdata())
    rest = np.squeeze(nb.load(input_rest).get_fdata())

    # mask data
    contrast_pos = contrast.copy()
    contrast_neg = contrast.copy()

    contrast_pos[contrast < min_contrast] = np.NaN
    contrast_neg[contrast > -min_contrast] = np.NaN

    contrast_pos[~label] = np.NaN
    contrast_neg[~label] = np.NaN
    contrast_pos[~np.isnan(contrast_pos)] = 1
    contrast_neg[~np.isnan(contrast_neg)] = 1

    # get resting-state data in mask
    rest_pos = rest.copy()
    rest_neg = rest.copy()

    rest_pos = rest_pos * contrast_pos
    rest_pos = rest_pos[~np.isnan(rest_pos)]
    rest_pos = rest_pos[rest_pos != np.min(rest_pos)]

    rest_neg = rest_neg * contrast_neg
    rest_neg = rest_neg[~np.isnan(rest_neg)]
    rest_neg = rest_neg[rest_neg != np.min(rest_neg)]

    # select random number of vertices
    label_shuffled_pos = random.sample(range(0, len(rest_pos)), nvert)
    label_shuffled_neg = random.sample(range(0, len(rest_neg)), nvert)

    rest_pos = rest_pos[label_shuffled_pos]
    rest_neg = rest_neg[label_shuffled_neg]

    # independent samples t-test
    # Levene's test is run to check for equal variances. If variances are not
    # equal, Welch's t-test is performed.
    _, p_levene = levene(rest_pos, rest_neg)
    if p_levene < 0.05:
        t, p = ttest_ind(rest_pos, rest_neg, equal_var=False)
    else:
        t, p = ttest_ind(rest_pos, rest_neg, equal_var=True)

    return rest_pos, rest_neg, t, p, p_levene


def analyze_fft(
    arr, fovx, fovy, xv, yv, f_cut=0.05, p_min=None, p_max=None, nsample=1000
):
    """This function computes the peak frequency and its corresponding power from a
    one-sided power spectrum sampled on a projection lines of a 2d array using nearest
    neighbor interpolation.

    Parameters
    ----------
    arr : ndarray
        2D input array.
    fovx : TYPE
        Field of view in x-direction (mm).
    fovy : float
        Field of view in y-firection (mm).
    xv : float
        x-coordinate of pca eigenvector.
    yv : float
        y-coordinate of pca eigenvector.
    f_cut : float, optional
        Cut off central spatial frequencies for peak detection. The default is
        0.05.
    p_min : float, optional
        Minimum prominence for peak detection. The default is None.
    p_max : float, optional
        Maximum prominence for peak detection. The default is None.
    nsample : int, optional
        Number of sampling point along projection line. The default is 1000.

    Returns
    -------
    k_max : float
        Spatial frequency with maximum peak power P_max.
    P_max : float
        Maximum peak power relative to central frequency along projection axis.
    k_line : ndarray
        Spatial frequencies along projection axis.
    fft_line : ndarray
        Corresponding spectral power along projection axis.

    """
    # add one if nsample is an odd integer
    if np.mod(nsample, 2) > 0:
        nsample += 1

    # get size of array
    x_size = np.shape(arr)[0]
    y_size = np.shape(arr)[1]

    # get coordinate of array
    kx = np.linspace(-1, 1, x_size) * (x_size / (2 * fovx))
    ky = np.linspace(-1, 1, y_size) * (y_size / (2 * fovy))
    ky_mesh, kx_mesh = np.meshgrid(ky, kx)

    # get projection line coordinates from pca eigenvector. Because we
    # interpolate using nearest neighbors the axis enc point is <size>-1.
    if np.abs(xv) < np.abs(yv):
        y_line = np.linspace(0, y_size - 1, nsample)
        x_line = xv / yv * y_line + y_size / 2
        x_line = x_line + x_size / 2 - x_line[int(nsample / 2)]
    else:
        x_line = np.linspace(0, x_size - 1, nsample)
        y_line = yv / xv * x_line + x_size / 2
        y_line = y_line + y_size / 2 - y_line[int(nsample / 2)]

    # check that no line steps over array border
    if np.any(y_line > y_size - 1):
        x_line[y_line > y_size - 1] = np.nan
        y_line[y_line > y_size - 1] = np.nan
        x_line = x_line[~np.isnan(x_line)]
        y_line = y_line[~np.isnan(y_line)]

    if np.any(y_line < 0):
        x_line[y_line < 0] = np.nan
        y_line[y_line < 0] = np.nan
        x_line = x_line[~np.isnan(x_line)]
        y_line = y_line[~np.isnan(y_line)]

    if np.any(x_line > x_size - 1):
        y_line[x_line > x_size - 1] = np.nan
        x_line[x_line > x_size - 1] = np.nan
        x_line = x_line[~np.isnan(x_line)]
        y_line = y_line[~np.isnan(y_line)]

    if np.any(y_line < 0):
        y_line[x_line < 0] = np.nan
        x_line[x_line < 0] = np.nan
        x_line = x_line[~np.isnan(x_line)]
        y_line = y_line[~np.isnan(y_line)]

    kxx_line = kx_mesh[np.round(x_line).astype(int), np.round(y_line).astype(int)]
    kyy_line = ky_mesh[np.round(x_line).astype(int), np.round(y_line).astype(int)]

    # get final k-axis
    k_line = np.sqrt(kxx_line**2 + kyy_line**2)

    # get fourier spectrum
    array_fft = get_fft(arr)
    fft_line = array_fft[np.round(x_line).astype(int), np.round(y_line).astype(int)]

    # get one-sided spectrum
    arg_null = int(np.argwhere(k_line == np.min(k_line))[-1])
    fft_line = fft_line[arg_null:]
    k_line = k_line[arg_null:]

    # normalize by central frequency
    fft_line = fft_line / fft_line[0] * 100

    # first central k-space line for peak detection
    fft_cut = fft_line.copy()
    k_cut = k_line.copy()

    fft_cut = fft_cut[k_cut > f_cut]
    k_cut = k_cut[k_cut > f_cut]

    # find peaks
    peak = find_peaks(fft_cut, prominence=(p_min, p_max))[0]

    if len(peak) < 1:
        p_max = np.nan
        k_max = np.nan
    else:
        p_max = fft_cut[peak[np.argwhere(fft_cut[peak] == np.max(fft_cut[peak]))[0]]]
        k_max = k_cut[peak[np.argwhere(fft_cut[peak] == np.max(fft_cut[peak]))[0]]]

    return k_max, p_max, k_line, fft_line


def get_pca(arr, fft_threshold=0.25):
    """This function computes the principal axes from the thresholded fourier spectrum
    of an input array. They are found by estimating the eigenvectors of the moments of
    inertia matrix I wich is defined by I = sum_i=0^N ( y_i^2 & -x_iy_i \\
    -x_iy_i & x_i^2 ). x_i and y_i are the coordinates of the remaining points of the
    thresholded Fourier spectrum. This procedure mainly follows [1].

    Parameters
    ----------
    arr : ndarray
        Input array.
    fft_threshold : float, optional
        Thresholded for getting the coordinates of the central part. The
        default is 0.25.

    Returns
    -------
    x_v1 : float
        x-coordinate of the major axis.
    y_v1 : float
        y-coordinate of the major axis.
    x_v2 : float
        x-coordinate of the minor axis.
    y_v2 : float
        y-coordinate of the minor axis.

    References
    -------
    .. [1] Borri, Marco, et al. A novel approach to evaluate spatial resolution
    of MRI clinical images for optimization and standardization of breast
    screening protocols, Med Phys 43(12), 6354--6363 (2016).

    """
    # compute normalized fourier spectrum of input array
    data_fft = get_fft(arr, write_output=False, normalization=True, n=10)

    # threshold spectrum
    data_fft[data_fft < fft_threshold] = 0
    data_fft[data_fft != 0] = 1

    # define mesh grid
    x = np.linspace(-np.shape(arr)[0] / 2, np.shape(arr)[0] / 2, np.shape(arr)[0])
    y = np.linspace(-np.shape(arr)[1] / 2, np.shape(arr)[1] / 2, np.shape(arr)[1])
    y_mesh, x_mesh = np.meshgrid(y, x)

    # get coordinates of thresholded spectrum
    x_mesh = x_mesh[data_fft != 0]
    y_mesh = y_mesh[data_fft != 0]

    # compute moments of inertia
    a = np.sum(x_mesh**2)
    b = np.sum(y_mesh**2)
    c = np.sum(x_mesh * y_mesh)

    m = np.array([[b, -c], [-c, a]])

    # compute eigenvalues and eigenvectors of matrix I
    evals, evecs = np.linalg.eig(m)

    # sort eigenvalues in decreasing order
    sort_indices = np.argsort(evals)[::-1]
    x_v1, y_v1 = evecs[:, sort_indices[0]]  # Eigenvector with largest eigenvalue
    x_v2, y_v2 = evecs[:, sort_indices[1]]

    return x_v1, y_v1, x_v2, y_v2


def get_fft(
    arr, write_output=False, path_output="", name_output="", normalization=False, n=1000
):
    """This function computes the power spectrum of a 2D numpy array. The result is
    saved as nifti image. Optionally, the FFT spectrum can be normalized. To do so, the
    input is shuffled N times and the 50th percentile of the power spectrum of the
    shuffled input array is computed. The mean 50th percentile of all shuffles is
    subtracted from the initial power spectrum and negative values are set to zero.
    Additionally, the resulting spectrum is normalized by its maximum power.

    Parameters
    ----------
    arr : ndarray
        2D nifti input array.
    write_output : bool, optional
        If output is written as nifti file. The default is False.
    path_output : str, optional
        Path where output is saved. The default is "".
    name_output : str, optional
        Basename of output image. The default is "".
    normalization : bool, optional
        Indicate if power spectrum is normalized. The default is False.
    n : int, optional
        Number of shuffled for normalization. The default is 1000.

    Returns
    -------
    array_fft : ndarray
        Fourier transformed array.

    """
    # copy input data
    data = arr.copy()

    # compute autocorrelation
    array_fft = np.abs(fftshift(fft2(data)))

    shuffle_mean = []
    if normalization is True:
        for i in range(n):
            # shuffle input array
            shuffle(data)  # shuffle along first axis
            data = data.T  # transpose array
            shuffle(data)  # shuffle along second axis
            data = data.T

            # get fft
            array_shuffle_fft = np.abs(fftshift(fft2(data)))

            # get 50th percentile
            shuffle_mean.append(np.percentile(array_shuffle_fft, 50))

        # get mean of percentiles
        array_50 = np.mean(shuffle_mean)

        # subtract 50th percentile from fft and divide by maximum power
        array_fft = array_fft - array_50
        array_fft[array_fft < 0] = 0
        array_fft = array_fft / np.max(array_fft)

    # write nifti
    if write_output is True:
        output = nb.Nifti1Image(array_fft, np.eye(4), nb.Nifti1Header())
        nb.save(output, os.path.join(path_output, name_output + "_fft.nii"))

    return array_fft


def get_bandpass(
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
    """This function generates a bandpass filter for an input image defined by lower and
    upper cutoff frequencies. Additionally, Gaussian attenuation of border can be appled
    and the filter can be only defined to a specific k-space region defined by lower and
    upper angle thresholds.

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


def get_rf_pulse_bw(file_in, npad=1000, ninterp=1000000):
    """This function computes the bandwidth in Hz of an RF pulse. The RF pulse shape has
    to be exported from the POET simulation. This can be done by clicking on the wanted
    pulse in the pulse sequence diagram and saving the event block as textfile. The
    bandwidth is calculated as the width at the threshold value of the normalized
    frequency magnitude.

    Parameters
    ----------
    file_in : str
        Textfile of one event block.
    npad : int, optional
        Pad zeros before and after pulse. The default is 1000.
    ninterp : int, optional
        Number of frequency steps for spline interpolation. The default is
        1000000.

    Returns
    -------
    rf : float
        Pulse shape in time domain (microseconds).
    rf_fft : float
        Pulse shape in frequency domain (cycles/second).
    bw : float
        Bandwidth (Hz).

    """
    # read file
    file = []
    with open(file_in) as fp:
        line = fp.readline()
        while line:
            line = fp.readline()
            file.append(line)

    # remove header from array
    file = file[9:]

    # get rf-pulse shape from array
    rf = []
    for i in range(len(file)):
        temp = [float(s) for s in re.findall(r"-?\d+\.?\d*", file[i])]
        if len(temp) > 0:
            rf.append(temp[0])

    # remove zeros at the beginning
    exit_loop = 0
    cnt = 0
    while exit_loop == 0:
        if rf[cnt] == 0:
            cnt += 1
        else:
            exit_loop = 1
    rf = rf[cnt:]

    # remove zeros at the end
    exit_loop = 0
    cnt = -1
    while exit_loop == 0:
        if rf[cnt] == 0:
            cnt -= 1
        else:
            exit_loop = 1
    rf = rf[: cnt + 1]

    # pad beginning and ending with zeros
    rf = np.pad(rf, npad, "constant")

    # get magnitude fft
    rf_fft = np.abs(fft(rf))

    # get frequency axis
    freq = np.fft.fftfreq(len(rf), 1e-6)  # timesteps in seconds

    # shift center frequency to center and normalize magnitude
    freq = fftshift(freq)
    rf_fft = fftshift(rf_fft)
    rf_fft = rf_fft / np.max(rf_fft)

    # get spline interpolation of fft
    freq_spline = np.linspace(np.min(freq), np.max(freq), ninterp)
    rf_fft_spline = np.interp(freq_spline, freq, rf_fft)

    threshold = 0.1
    for i in range(len(rf_fft_spline) - 1):
        if rf_fft_spline[i] <= threshold < rf_fft_spline[i + 1]:
            freq1 = freq_spline[i]

        if rf_fft_spline[i] > threshold >= rf_fft_spline[i + 1]:
            freq2 = freq_spline[i]

    # bandwidth in Hz
    bw = freq2 - freq1
    print("Pulse bandwidth in Hz: " + str(bw))

    # plot pulse and corresponding fft
    fig, ax = plt.subplots()
    ax.plot(np.linspace(0, len(rf), len(rf)), rf)
    ax.set_xlabel("Time in microseconds")
    ax.set_ylabel("Amplitude in a.u.")
    ax.set_title("Pulse shape in time domain")
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(freq_spline, rf_fft_spline)
    ax.set_xlabel("frequency in cycles/second")
    ax.set_ylabel("Normalized amplitude in a.u.")
    ax.set_title("Pulse shape in frequency domain")
    plt.show()

    return rf, rf_fft, bw


def get_bandwidth_per_px(file_in, file_out):
    """This function reads the bandwidth per pixel phase encode from siemens dicom data.
    This parameter is needed for fieldmap undistortion.

    Parameters
    ----------
    file_in : str
        File name of dicom file.
    file_out : str
        File name of output text file.

    Raises
    ------
    ValueError
        If input file does not have the right file extension.
    """
    # make output folder
    out = get_filename(file_out)
    if not os.path.exists(out[0]):
        os.makedirs(out[0])

    # read dicom
    with open(file_out, "w", encoding="utf-8") as f:
        if not (file_in.endswith(".ima") or file_in.endswith(".dcm")):
            raise ValueError("Input file should have *.ima or *.dcm extension!")
        data = pydicom.dcmread(file_in)
        if ["0019", "1028"] in data:
            f.write(
                f"{data.SeriesNumber}\t{data.SeriesDescription}\t{data['0019', '1028']}"
            )
        else:
            print("Dicom tag not found.")
