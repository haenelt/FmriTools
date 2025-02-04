# -*- coding: utf-8 -*-
"""Time-series manipulation."""

import os
from math import prod

import nibabel as nb
import numpy as np
from nilearn.signal import clean
from numpy.fft import fft, ifft
from scipy.interpolate import InterpolatedUnivariateSpline as Interp
from scipy.ndimage import gaussian_filter

from ..io.filename import get_filename
from .. import execute_command

__all__ = [
    "ScaleTimeseries",
    "FilterTimeseries",
    "slice_timing_correction",
    "interpolate",
    "average",
    "bandpass_afni",
]


class ScaleTimeseries:
    """Implementation of several methods to scale fmri time series data.

    Parameters
    ----------
    arr : np.ndarray
        4D time series array.

    """

    def __init__(self, arr):
        self.arr = arr
        self.nt = np.shape(arr)[3]  # number of time points

    @property
    def arr_mean(self):
        """Temporal mean."""
        return np.nanmean(self.arr, axis=3)

    @property
    def arr_std(self):
        """Temporal standard deviation."""
        return np.nanstd(self.arr, axis=3)

    @property
    def arr_mean_repeated(self):
        """Expanded temporal mean."""
        return np.repeat(self.arr_mean[:, :, :, np.newaxis], self.nt, axis=3)

    @property
    def arr_std_repeated(self):
        """Expanded temporal standard deviation."""
        return np.repeat(self.arr_std[:, :, :, np.newaxis], self.nt, axis=3)

    def psc(self, cutoff_size=None):
        """Percent signal change conversion."""
        self.arr = self._save_division(self.arr, self.arr_mean_repeated) * 100
        if cutoff_size:
            self.arr = self._cutoff(100, cutoff_size)
        return self.arr

    def normalize(self, cutoff_size=None):
        """Normalize time series."""
        self.arr = self._save_division(self.arr, self.arr_mean_repeated)
        if cutoff_size:
            self.arr = self._cutoff(1, cutoff_size)
        return self.arr

    def standardize(self, cutoff_size=None):
        """Standardize time series."""
        self.arr = self._save_division(
            self.arr - self.arr_mean_repeated, self.arr_std_repeated
        )
        if cutoff_size:
            self.arr = self._cutoff(0, cutoff_size)
        return self.arr

    def demean(self, cutoff_size=None):
        """Demean time series."""
        self.arr = self._save_division(
            self.arr - self.arr_mean_repeated, self.arr_mean_repeated
        )
        if cutoff_size:
            self.arr = self._cutoff(0, cutoff_size)
        return self.arr

    def _cutoff(self, mu, width):
        """Threshold outliers."""
        arr_max = np.max(self.arr, axis=3)
        arr_min = np.min(self.arr, axis=3)
        self.arr[arr_max > mu + width, :] = mu + width
        self.arr[arr_min < mu - width, :] = mu - width
        return self.arr

    @staticmethod
    def _save_division(arr1, arr2):
        return np.divide(arr1, arr2, out=np.zeros_like(arr1), where=arr2 != 0)

    @classmethod
    def from_file(cls, file_data):
        """Initialize class object from file."""
        data = nb.load(file_data)
        return cls(data.get_fdata())


class FilterTimeseries(ScaleTimeseries):
    """Implementation of different filters to low- or highpass fmri time series data.

    Parameters
    ----------
    arr : np.ndarray
        4D time series array.
    TR : float
        Repetition time in s.

    """

    def __init__(self, arr, TR):
        super().__init__(arr)
        self.TR = TR

    def detrend(self, cutoff_sec, store_dc=False):
        """Detrend time series by convolving the time series with a gaussian running
        line smoother [1]_. This is also the method used by FSL for detrending.

        Returns
        -------
        cutoff_sec : float
            Highpass 1/cutoff frequency in Hz.
        store_dc : bool
            Add dc component.

        Returns
        -------
        ndarray
            Filtered array.

        References
        -------
        .. [1] https://lukas-snoek.com/NI-edu/fMRI-introduction/week_4/temporal_preprocessing.html

        """
        arr0 = self.arr_mean_repeated  # dc component
        sigma = cutoff_sec / (np.sqrt(8 * np.log(2)) * self.TR)  # kernel definition
        arr_filtered = np.zeros_like(self.arr)
        for ix, iy, iz in np.ndindex(self.arr.shape[:-1]):
            arr_filtered[ix, iy, iz] = gaussian_filter(self.arr[ix, iy, iz], sigma)
        self.arr = self.arr - arr_filtered  # remove lowpass
        # store dc component
        if store_dc:
            self.arr = self.arr + arr0
        return self.arr

    def lowpass_sma(self, window_size):
        """Filters time series data using a simple moving average (SMA) filter. The SMA
        filter is a simple average of the data points within a user-defined window. The
        window is centered on each data point. Note that an even filter of size n will
        effectively have the same filter size as an odd filter of size n+1.

        Parameters
        ----------
        window_size : int
            Size of the window.

        Returns
        -------
        ndarray
            Filtered array.

        """
        # parameters
        n = int(window_size / 2)  # one-sided window size

        arr2d = self._reshape(self.arr)
        arr2d_filt = np.zeros_like(arr2d)
        for t in range(self.nt):
            if t < n:
                arr2d_filt[:, t] = np.mean(arr2d[:, : t + n], axis=1)
            elif t > self.arr.shape[-1] - n:
                arr2d_filt[:, t] = np.mean(arr2d[:, t - n :], axis=1)
            else:
                arr2d_filt[:, t] = np.mean(arr2d[:, t - n : t + n], axis=1)

        return self._reshape_back(arr2d_filt, np.shape(self.arr))

    def lowpass_gaussian(self, sigma, normalize=False):
        """Filters time series data using a Gaussian filter. Thr filter size can be
        defined by the standard deviation of the Gaussian kernel.

        Parameters
        ----------
        sigma : float
            Standard deviation of the Gaussian function in s.
        normalize : bool, optional
            Normalize the filter to max(filter) = 1.

        Returns
        -------
        ndarray
            Filtered array.

        """
        # gaussian filter definition
        x = np.linspace(-int(self.nt / 2), int(self.nt / 2), self.nt) * self.TR
        f = [
            1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-(t**2) / (2 * sigma**2))
            for t in x
        ]
        f_fft = np.abs(fft(f))

        # normalize filter
        if normalize:
            f_fft /= np.max(f_fft)

        # apply filter in spatial frequency space
        arr2d_filt = self._apply_filter(self.arr, f_fft)

        return self._reshape_back(arr2d_filt, np.shape(self.arr))

    def bandpass_boxcar(self, cutoff_low=10, cutoff_high=1000, preserve_range=False):
        """Filters time series data using a boxcar filter. Optionally, only the lowpass
        or highpass filter can be applied by setting the other cutoff to None. The code
        is adapted from [1]_. Optionally, the original range of the data can be
        preserved by keeping the spatial center frequency in the filter.

        Parameters
        ----------
        cutoff_low : float, optional
            Cutoff frequency for the lowpass filter in 1/Hz.
        cutoff_high : float, optional
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
        fs = 1.0 / self.TR

        lowidx = int(self.nt / 2) + 1
        if cutoff_low is not None:
            lowpass_freq = 1.0 / cutoff_low
            lowidx = np.round(lowpass_freq / fs * self.nt).astype(int)

        highidx = 0
        if cutoff_high is not None:
            highpass_freq = 1.0 / cutoff_high
            highidx = np.round(highpass_freq / fs * self.nt).astype(int)

        f_fft = np.zeros(self.nt)
        f_fft[highidx:lowidx] = 1
        if preserve_range:
            f_fft[0] = 1.0

        f_fft = ((f_fft + f_fft[::-1]) > 0).astype(int)
        if np.all(f_fft == 1):
            print("No filtering applied.")
            return self.arr

        # apply filter in spatial frequency space
        arr2d_filt = self._apply_filter(self.arr, f_fft)

        return self._reshape_back(arr2d_filt, np.shape(self.arr))

    def bandpass_butterworth(self, cutoff_low=10, cutoff_high=1000):
        """Filters time series data using a Butterworth filter. Optionally, only the
        lowpass or highpass filter can be applied by setting the other cutoff to None.

        Parameters
        ----------
        cutoff_low : float, optional
            Cutoff frequency for the lowpass filter in 1/Hz.
        cutoff_high : float, optional
            Cutoff frequency for the highpass filter in 1/Hz.

        Returns
        -------
        ndarray
            Filtered array.

        """
        # parameters
        low_pass = 1.0 / cutoff_low if cutoff_low is not None else None
        high_pass = 1.0 / cutoff_high if cutoff_high is not None else None

        arr2d = self._reshape(self.arr)
        arr2d_filt = np.zeros_like(arr2d)
        ndata = arr2d.shape[0]
        for i in range(ndata):
            # filter parameters can be butterworth, cosine (lowpass not available)
            # or False. Optionally, the time series can be standardized (zscore) or
            # converted to a percent signal change (psc).
            arr2d_filt[i, :] = clean(
                arr2d[i, :],
                filter="butterworth",
                detrend=False,
                standardize=False,
                t_r=self.TR,
                low_pass=low_pass,
                high_pass=high_pass,
            )

        return self._reshape_back(arr2d_filt, np.shape(self.arr))

    def _apply_filter(self, arr, f_fft):
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
        arr2d = self._reshape(arr)
        arr2d_filt = np.zeros_like(arr2d)
        ndata = arr2d.shape[0]
        for i in range(ndata):
            arr2d_filt[i, :] = np.real(ifft(fft(arr2d[i, :]) * f_fft))

        return arr2d_filt

    @staticmethod
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

    @staticmethod
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


def slice_timing_correction(
    file_in, TR_old, TR_new, order, mb=None, manufacturer="siemens", prefix="a"
):
    """This function performs slice timing correction of a nifti time series. For
    interleaved slice ordering, interleaved ascending is assumed. The correction is done
    by temporal interpolation of single voxel time series using cubic interpolation. To
    omit extrapolation errors at the edges, the first and last volumes of the time
    series are appended at the beginning and at the end, respectively. These time points
    are removed again after the interpolation step. The interpolated time series is
    sampled onto a regular grid with a defined new TR. Therefore, the reference slice is
    always the first slice acquired at t = 0. For time series acquired with multiband,
    the number of slices has to be a multiple of the multiband factor. For interleaved
    slice acquisition, siemens sequences start with the odd (even) slice for images with
    odd (even) number of slices. You can switch to cmrr to consider that the cmrr
    sequence always starts with the odd slice (first slice) in interleaved mode
    irrespective of the number of slices.

    Parameters
    ----------
    file_in : str
        Filename of nifti time series.
    TR_old : float
        TR of time series in seconds.
    TR_new : float
        TR of slice timing corrected time series in s.
    order : str
        Slice ordering (ascending, descending, interleaved).
    mb : int
        Multiband factor.
    manufacturer : str
        Sequence type (siemens, cmrr).
    prefix : str, optional
        Prefix of output time series basename. The default is "a".
    """
    if not mb:
        mb = 1

    if manufacturer not in ["siemens", "cmrr"]:
        raise ValueError("Unknown manufacturer!")

    # get filename
    path_file, name_file, ext_file = get_filename(file_in)

    # load data
    data = nb.load(file_in)
    nz = data.header["dim"][3]

    # effective number of sequentially acquired slices
    mb_package = nz / mb
    if np.mod(mb_package, 1):
        raise ValueError("Number of slices and multiband factor does not match!")
    else:
        mb_package = int(mb_package)

    # spatial order of acquired slices
    if order == "ascending":
        slice_order = np.arange(0, nz)
    elif order == "descending":
        slice_order = np.arange(nz - 1, -1, -1)
    elif (
        order == "interleaved" and np.mod(nz, 2) or manufacturer == "cmrr"
    ):  # odd slice number
        slice_order = np.arange(0, nz, 2)
        slice_order = np.append(slice_order, np.arange(1, nz, 2))
    elif order == "interleaved" and not np.mod(nz, 2):  # even slice number
        slice_order = np.arange(1, nz, 2)
        slice_order = np.append(slice_order, np.arange(0, nz, 2))
    else:
        raise ValueError("Choose a valid slice ordering!")

    # temporal order of acquired slices
    if order == "interleaved":
        target = np.ceil(mb_package / 2).astype(int)
        temporal_order = np.arange(0, target)
        if float(mb_package / 2).is_integer():
            temporal_order2 = np.arange(0, target)
        else:
            temporal_order2 = np.arange(0, target)[:-1]

        temporal_order = np.tile(temporal_order, mb)
        temporal_order2 = np.tile(temporal_order2, mb)
        temporal_order = np.append(temporal_order, temporal_order2 + target)
    else:
        temporal_order = np.arange(0, mb_package)
        temporal_order = np.tile(temporal_order, mb)

    # some prints for sanity check
    print("Spatial order of slices: " + str(slice_order))
    print("Temporal order of slices: " + str(temporal_order))

    # temporal interpolation
    TA = TR_old / mb_package  # acquisition time needed for one slice
    sequence = np.column_stack((slice_order, temporal_order))
    file_out = os.path.join(path_file, prefix + name_file + ext_file)
    interpolate(file_in, file_out, TR_old, TR_new, TA, sequence)


def interpolate(file_in, file_out, tr_old, tr_new, ta=None, sequence=None):
    """This function performs temporal interpolation of voxel time courses using cubic
    interpolation and samples data with a new repetition time.

    Parameters
    ----------
    file_in : str
        Filename of nifti time series.
    file_out : str
        Filename of output nifti time series.
    tr_old : float
        TR of time series in seconds.
    tr_new : float
        TR of interpolated time series in s.
    ta : float
        Time for a single slice (z-direction) to acknowledge multiband acquisitions.
        ta=tr_old is not set. The default is "None".
    sequence : str
        Sequence of acquired slices. This parameter is a numpy array with dimensions
        (number of slices, 2). sequence(:, 0) indicates the sampled z-slice and
        sequence(:, 1) the corresponding time. This parameter acknowledges differences
        in slice timinings in 2D acqisitions. If not set, it is assumed that all slices
        were acquired at the same time. The default is "None".
    """
    # load data and temporarilly add volumes at the beginning and add for cubic
    # interpolation
    data = nb.load(file_in)
    nx, ny, nz, nt = data.header["dim"][1:5]
    arr = np.zeros((nx, ny, nz, nt + 2))
    arr[:, :, :, 0] = data.get_fdata()[:, :, :, 0]
    arr[:, :, :, -1] = data.get_fdata()[:, :, :, -1]
    arr[:, :, :, 1:-1] = data.get_fdata()

    # default parameters
    if ta is None:
        ta = tr_old
    if sequence is None:
        sequence = np.ones((nz, 2), dtype=int)
        sequence[:, 0] = np.arange(nz)

    tt = tr_old * nt  # total acquisition time
    tr_append = np.floor(tr_old / tr_new).astype(int) * tr_new  # number of appended TRs
    t_new = np.arange(-tr_append, tt + tr_append, tr_new)  # grid points of output array

    # loop over voxels and apply cubic interpolation
    arr_corrected = np.zeros((nx, ny, nz, len(t_new)))
    for z in range(nz):
        for x in range(nx):
            for y in range(ny):
                t = np.arange(
                    sequence[z, 1] * ta - tr_old,
                    sequence[z, 1] * ta + (nt + 1) * tr_old,
                    tr_old,
                )
                t = t[: nt + 2]  # ensure same size of arrays
                cubic_interper = Interp(t, arr[x, y, int(sequence[z, 0]), :], k=3)
                arr_corrected[x, y, int(sequence[z, 0]), :] = cubic_interper(t_new)

    # delete appended volumes
    vols_keep1 = t_new >= 0
    vols_keep2 = t_new < tt
    vols_keep = vols_keep1 * vols_keep2
    arr_corrected = arr_corrected[:, :, :, vols_keep]

    # clean corrected data
    arr_min = np.min(data.get_fdata())
    arr_max = np.max(data.get_fdata())
    arr_corrected[np.isnan(arr_corrected)] = 0.0
    arr_corrected[arr_corrected < arr_min] = 0.0
    arr_corrected[arr_corrected > arr_max] = arr_max

    # update header
    data.header["dim"][4] = np.shape(arr_corrected)[3]
    data.header["datatype"] = 16

    # write output
    output = nb.Nifti1Image(arr_corrected, data.affine, data.header)
    output = _set_tr(output, tr_new)
    nb.save(output, file_out)


def _set_tr(img, tr):
    """Helper function to set tr in nifti header."""
    header = img.header.copy()
    zooms = header.get_zooms()[:3] + (tr,)
    header.set_zooms(zooms)

    return img.__class__(img.get_fdata().copy(), img.affine, header)


def average(img_input, path_output, name_output):
    """This function computes the element-wise average of multiple nifti files.

    Parameters
    ----------
    img_input : str
        List of nifti input paths.
    path_output : str
        Path where output is saved.
    name_output : str
        Basename of output files.
    """
    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # load first dataset to initialize final time series
    res = nb.load(img_input[0])
    res_array = np.zeros_like(res.get_fdata())

    for i, _ in enumerate(img_input):
        res_array += nb.load(img_input[i]).get_fdata()

    # divide summed time series by number of time series
    res_array = res_array / len(img_input)

    # write output
    output = nb.Nifti1Image(res_array, res.affine, res.header)
    nb.save(output, os.path.join(path_output, name_output + ".nii"))


def bandpass_afni(file_in, file_out, TR, lp_freq, hp_freq):
    """Bandpass filter time series using AFNI.

    Parameters
    ----------
    file_in : str
        File name of input file.
    file_out : str
        File name of output file.
    TR : float
        Repetition time in s.
    lp_freq : float
        Lowpass filter frequencyin Hz.
    hp_freq : float
        Highpass filter frequency in Hz.
    """
    command = "3dBandpass"
    command += f" -prefix {file_out}"
    command += f" -dt {TR:.6f}"
    command += f" {hp_freq:.6f} {lp_freq:.6f}"
    command += f" {file_in}"

    # run
    execute_command(command)
