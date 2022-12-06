# -*- coding: utf-8 -*-

# python standard library inputs
from math import prod

# external inputs
import numpy as np
import nibabel as nb
from numpy.fft import fft, ifft
from scipy.ndimage import gaussian_filter
from nilearn.signal import clean

__all__ = ["ScaleTimeseries", "FilterTimeseries"]


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
                arr2d_filt[:, t] = np.mean(arr2d[:, t - n:], axis=1)
            else:
                arr2d_filt[:, t] = np.mean(arr2d[:, t - n: t + n], axis=1)

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
