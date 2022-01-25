# -*- coding: utf-8 -*-

# python standard library inputs
from math import prod

# external inputs
import numpy as np
from numpy.fft import fft, ifft
from nilearn.signal import clean

__all__ = ['lowpass_sma', 'lowpass_gaussian', 'bandpass_boxcar',
           'bandpass_butterworth']


def lowpass_sma(arr, window_size):
    """Running mean filter.

    Filters time series data using a simple moving average (SMA) filter. The SMA
    filter is a simple average of the data points within a user-defined window.
    The window is centered on each data point. Note that an even filter of size
    n will effectively have the same filter size as an odd filter of size n+1.

    Parameters
    ----------
    arr : ndarray
        Array to be filtered. Time must be the last dimension.
    window_size : int
        Size of the window.

    Returns
    -------
    ndarray
        Filtered array.

    """

    # parameters
    n = int(window_size / 2)  # one-sided window size
    nt = arr.shape[-1]  # number of time points

    arr2d = _reshape(arr)
    arr2d_filt = np.zeros_like(arr2d)
    for t in range(nt):
        if t < n:
            arr2d_filt[:, t] = np.mean(arr2d[:, :t+n], axis=1)
        elif t > arr.shape[-1] - n:
            arr2d_filt[:, t] = np.mean(arr2d[:, t-n:], axis=1)
        else:
            arr2d_filt[:, t] = np.mean(arr2d[:, t-n:t+n], axis=1)

    return _reshape_back(arr2d_filt, np.shape(arr))


def lowpass_gaussian(arr, tr, sigma, normalize=False):
    """Gaussian lowpass filter.

    Filters time series data using a Gaussian filter. Thr filter size can be
    defined by the standard deviation of the Gaussian kernel.

    Parameters
    ----------
    arr : ndarray
        Array to be filtered. Time must be the last dimension.
    tr : float
        Repetition time in s.
    sigma : float
        Standard deviation of the Gaussian function in s.
    normalize : bool, optional
        Normalize the filter to max(filter) = 1.

    Returns
    -------
    ndarray
        Filtered array.

    """

    # parameters
    nt = arr.shape[-1]

    # gaussian filter definition
    filter_range = np.linspace(-int(nt / 2), int(nt / 2), nt) * tr
    f = [
        1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-t ** 2 / (2 * sigma ** 2))
        for t in filter_range]
    f_fft = np.abs(fft(f))

    # normalize filter
    if normalize:
        f_fft /= np.max(f_fft)

    # apply filter in spatial frequency space
    arr2d_filt = _apply_filter(arr, f_fft)

    return _reshape_back(arr2d_filt, np.shape(arr))


def bandpass_boxcar(arr, tr, cutoff_lowpass=10, cutoff_highpass=1000,
                    preserve_range=False):
    """Boxcar bandpass filter.

    Filters time series data using a boxcar filter. Optionally, only the lowpass
    or highpass filter can be applied by setting the other cutoff to None. The
    code is adapted from [1]_. Optionally, the original range of the data can be
    preserved by keeping the spatial center frequency in the filter.

    Parameters
    ----------
    arr : ndarray
        Array to be filtered. Time must be the last dimension.
    tr : float
        Repetition time in s.
    cutoff_lowpass : float, optional
        Cutoff frequency for the lowpass filter in 1/Hz.
    cutoff_highpass : float, optional
        Cutoff frequency for the highpass filter in 1/Hz.
    preserve_range : bool, optional
        Preserve the range of the data.

    Returns
    -------
    ndarray
        Filtered array.

    References
    ----------
    .. [1] https://nipype.readthedocs.io/en/latest/users/examples/rsfmri_vol_surface_preprocessing.html

    """

    # parameters
    fs = 1.0 / tr
    nt = arr.shape[-1]

    lowidx = int(nt / 2) + 1
    if cutoff_lowpass is not None:
        lowpass_freq = 1.0 / cutoff_lowpass
        lowidx = np.round(lowpass_freq / fs * nt).astype(int)

    highidx = 0
    if cutoff_highpass is not None:
        highpass_freq = 1.0 / cutoff_highpass
        highidx = np.round(highpass_freq / fs * nt).astype(int)

    f_fft = np.zeros(nt)
    f_fft[highidx:lowidx] = 1

    if preserve_range:
        f_fft[0] = 1.0

    f_fft = ((f_fft + f_fft[::-1]) > 0).astype(int)

    if np.all(f_fft == 1):
        print("No filtering applied.")
        return arr

    # apply filter in spatial frequency space
    arr2d_filt = _apply_filter(arr, f_fft)

    return _reshape_back(arr2d_filt, np.shape(arr))


def bandpass_butterworth(arr, tr, cutoff_lowpass=10, cutoff_highpass=1000):
    """Butterworth bandpass filter.

    Filters time series data using a Butterworth filter. Optionally, only the
    lowpass or highpass filter can be applied by setting the other cutoff to
    None.

    Parameters
    ----------
    arr : ndarray
        Array to be filtered. Time must be the last dimension.
    tr : float
        Repetition time in s.
    cutoff_lowpass : float, optional
        Cutoff frequency for the lowpass filter in 1/Hz.
    cutoff_highpass : float, optional
        Cutoff frequency for the highpass filter in 1/Hz.

    Returns
    -------
    ndarray
        Filtered array.

    """

    # parameters
    low_pass = 1.0 / cutoff_lowpass if cutoff_lowpass is not None else None
    high_pass = 1.0 / cutoff_highpass if cutoff_highpass is not None else None

    arr2d = _reshape(arr)
    arr2d_filt = np.zeros_like(arr2d)
    ndata = arr2d.shape[0]
    for i in range(ndata):
        # filter parameters can be butterworth, cosine (lowpass not available)
        # or False. Optionally, the time series can be standardized (zscore) or
        # converted to a percent signal change (psc).
        arr2d_filt[i, :] = clean(arr2d[i, :],
                                 filter="butterworth",
                                 detrend=False,
                                 standardize=False,
                                 t_r=tr,
                                 low_pass=low_pass,
                                 high_pass=high_pass
                                 )

    return _reshape_back(arr2d_filt, np.shape(arr))


def _apply_filter(arr, f_fft):
    """Apply filter in spatial frequency space.

    Parameters
    ----------
    arr : ndarray
        Array to be filtered.
    f_fft : ndarray
        Filter in spatial frequency space.

    Returns
    -------
    ndarray
        Filtered array.

    """

    arr2d = _reshape(arr)
    arr2d_filt = np.zeros_like(arr2d)
    ndata = arr2d.shape[0]
    for i in range(ndata):
        arr2d_filt[i, :] = np.real(ifft(fft(arr2d[i, :]) * f_fft))

    return arr2d_filt


def _reshape(array):
    """Reshape array to 2D.

    Parameters
    ----------
    array : ndarray
        Array to be reshaped.

    Returns
    -------
    ndarray
        Array reshaped to 2D.

    """

    dim1 = prod(np.shape(array)[:-1])
    dim2 = np.shape(array)[-1]

    return np.reshape(array, (dim1, dim2))


def _reshape_back(array, shape):
    """Reshape array back to original shape.

    Parameters
    ----------
    array : ndarray
        Array to be reshaped.
    shape : tuple
        Shape of the original array.

    Returns
    -------
    ndarray
        Array reshaped to original shape.

    """

    return np.reshape(array, shape)
