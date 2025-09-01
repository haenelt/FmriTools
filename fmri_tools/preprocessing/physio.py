# -*- coding: utf-8 -*-
"""Utilities for time series noise removal."""

import copy
from pathlib import Path

import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
from nilearn.maskers import NiftiMasker
from scipy.ndimage import binary_erosion
from sklearn.decomposition import PCA

__all__ = ["CleanTimeseries", "fd_power", "friston_24", "global_signal", "compcor"]


class CleanTimeseries:
    """Estimate noise regressors and clean fMRI time series by regressing out confounds.
    Optionally detrend, smooth, and bandpass filter time series.

    Noise regressors can be estaimted based on:
        1) motion (Friston 24),
        2) CompCor,
        3) scrubbing,
        4) global signal.

    Parameters
    ----------
    fmri_img : nb.nifti1.Nifti1Image
        fMRI time series image.
    brain_img : nb.nifti1.Nifti1Image
        Brain mask image.
    tr : float
        Repetition time in s.
    n_erode : int, optional
        Erosion iterations for masks. Defaults to 3.

    """

    def __init__(self, fmri_img, brain_mask, tr: float, n_erode: int = 3):
        self.fmri_img = fmri_img
        self.n_erode = n_erode
        self.brain_mask = self._erode(brain_mask)
        self.tr = tr
        # Initialize confound array where all regressors will be stored.
        self.confounds = np.empty((fmri_img.shape[-1], 0))  # (n_timepoints, 0 columns)

    def __call__(
        self,
        motion,
        wm_mask,
        csf_mask=None,
        n_compcor=5,
        fd_threshold=0.5,
        use_global=False,
        high_pass=0.009,
        low_pass=0.08,
        fwhm=None,
        detrend=False,
    ):
        """Estimate confounds and clean fMRI data.

        Parameters
        ----------
        motion : np.ndarray
            Motion parameters of shape (T, 6).
        wm_mask : nb.nifti1.Nifti1Image
            Nibabel object containing a corresponding white matter mask.
        csf_mask : nb.nifti1.Nifti1Image, optional
            Nibabel object containing a corresponding CSF mask.
        n_compcor : int, optional
            Number of PCA components per tissue class.
        fd_threshold : float, optional
            FD power threshold in mm.
        use_global : bool, optional
            Add global signal confound.
        high_pass : float | None
            Temporal highpass filter in Hz (no filtering if set to None). Defaults to 
            0.009.
        low_pass : float | None
            Temporal lowpass filter in Hz (no filtering is set to None). Defaults to 
            0.08.
        fwhm : float | None
            Spatial gaussian filter size (in mm) for smoothing cleaned fMRI data.
        detrend : bool | None
            Apply detrending to cleaned fMRI data.

        """
        obj_copy = copy.deepcopy(self)
        obj_copy.add_motion(motion)
        obj_copy.add_compcor(wm_mask, csf_mask, n_compcor=n_compcor)
        obj_copy.add_scrubbing(motion, fd_threshold=fd_threshold)
        if use_global:
            obj_copy.add_global()

        return obj_copy.fit(high_pass, low_pass, fwhm, detrend)

    def add_motion(self, motion):
        """Add motion parameters and derived regressors in the Friston style to the
        list of confounds.

        Parameters
        ----------
        see `fmri_tools.preprocessing.physio.friston_24`

        """
        self.confounds = np.hstack((self.confounds, friston_24(motion)))

    def add_compcor(self, wm_mask, csf_mask=None, n_compcor=5):
        """Add compcor regressors to the list of confounds.

        Parameters
        ----------
        see `fmri_tools.preprocessing.physio.compcor`

        """
        _wm_mask = self._erode(wm_mask)
        _csf_mask = self._erode(csf_mask) if csf_mask is not None else None
        self.confounds = np.hstack(
            (self.confounds, compcor(self.fmri_img, _wm_mask, _csf_mask, n_compcor))
        )

    def add_scrubbing(
        self, motion, fd_threshold=0.5, rot_to_mm=50, rotations_in_degrees=True, max_gap=5
    ):
        """Add one-hot scrubbing regressors based on FD power estimates. Binary
        scrubbing regressors are computed by tresholding FD power traces with
        `fd_threshold` (in mm). If two outliers are separated by <= `max_gap` frames,
        the in-between time points are also flagged as outliers.

        Parameters
        ----------
        see `fmri_tools.preprocessing.physio.fd_power`

        """
        fd = fd_power(motion, rot_to_mm, rotations_in_degrees)
        n_tp = len(fd)

        # Initial outlier indices
        bad_idx = np.where(fd > fd_threshold)[0]
        if len(bad_idx) == 0:
            return
        
        # Expand to include time points between close outliers
        expanded_idx = set()
        prev = bad_idx[0]
        expanded_idx.add(prev)

        for idx in bad_idx[1:]:
            if idx - prev <= max_gap:
                # include all intermediate points
                expanded_idx.update(range(prev, idx + 1))
            expanded_idx.add(idx)
            prev = idx

        expanded_idx = sorted(expanded_idx)

        # One-hot regressors
        scrub = np.zeros((n_tp, len(expanded_idx)), dtype=int)
        scrub[expanded_idx, np.arange(len(expanded_idx))] = 1

        # Append to confounds
        self.confounds = np.hstack((self.confounds, scrub))

    def add_global(self):
        """Add global signal regressors to the list of confounds."""
        self.confounds = np.hstack(
            (self.confounds, global_signal(self.fmri_img, self.brain_mask))
        )

    def fit(
        self,
        standardize=False,
        high_pass=0.009,
        low_pass=0.08,
        fwhm=None,
        detrend=False,
    ):
        """Clean fmri be regressing out confound design matrix. Optionally, cleaned
        times series will be detrended, spatially smoothed and bandpass filtered.

        Default bandpass filtering parameters are taken from Power et al., Neuroimage
        2014, Laumann et al., Neuron 2015.

        Parameters
        ----------
        standardize : bool
            Standardize cleaned fMRI time series.
        high_pass : float | None
            Temporal highpass filter in Hz (no filtering if set to None). Defaults to 
            0.009.
        low_pass : float | None
            Temporal lowpass filter in Hz (no filtering is set to None). Defaults to 
            0.08.
        fwhm : float | None
            Spatial gaussian filter size (in mm) for smoothing cleaned fMRI data.
        detrend : bool | None
            Apply detrending to cleaned fMRI data.

        """
        # Check if confounding design matrix was estimated.
        if len(self.confounds) == 0:
            raise ValueError("No confounds could be found!")

        masker = NiftiMasker(
            mask_img=self.brain_mask,
            standardize=standardize,
            standardize_confounds=True,
            t_r=self.tr,
            high_pass=high_pass,
            low_pass=low_pass,
            smoothing_fwhm=fwhm,
            detrend=detrend,
        )
        cleaned = masker.fit_transform(self.fmri_img, confounds=self.confounds)
        return masker.inverse_transform(cleaned)

    def plot_confound(self, fname):
        """Save plot of confound matrix.

        Parameters
        ----------
        fname : str
            File name of image saved to disk.

        """
        # Check if confounding design matrix was estimated.
        if len(self.confounds) == 0:
            raise ValueError("No confounds could be found!")

        dir_out = Path(fname).parent
        dir_out.mkdir(exist_ok=True, parents=True)

        # Standardize confound matrix
        _confounds = self.confounds.copy()
        _mean = np.mean(_confounds, axis=0, keepdims=True)
        _std = np.std(_confounds, axis=0, keepdims=True)
        _confounds = (_confounds - _mean) / _std

        plt.figure(figsize=(12, 6))

        # Plot heatmap
        im = plt.imshow(
            _confounds,
            aspect="auto",
            cmap="coolwarm",
            interpolation="nearest",
            vmin=-np.max(np.abs(self.confounds)),
            vmax=np.max(np.abs(self.confounds)),
        )

        plt.colorbar(im, label="Amplitude")
        plt.xlabel("Regressors")
        plt.ylabel("Timepoints")
        plt.title("Confound Design Matrix")
        plt.tight_layout()

        # Save figure
        plt.savefig(fname, dpi=300, bbox_inches="tight")

    def _erode(self, data, n_erode=None):
        """Apply binary erosion to nifti object."""
        n_erode = self.n_erode if n_erode is None else n_erode
        arr = data.get_fdata()
        arr_eroded = binary_erosion(arr, iterations=self.n_erode)
        return nb.Nifti1Image(arr_eroded.astype(np.uint8), data.affine, data.header)


def compcor(fmri_img, wm_mask, csf_mask=None, n_compcor=5):
    """Generate regressor from white matter and csf masks based on the CompCor approach.
    Only uses white matter mask if no csf mask is provided.

    Parameters
    ----------
    fmri_img : nb.nifti1.Nifti1Image
        Nibabel object containing an fMRI time series.
    wm_mask : nb.nifti1.Nifti1Image
        Nibabel object containing a corresponding white matter mask.
    csf_mask : nb.nifti1.Nifti1Image, optional
        Nibabel object containing a corresponding CSF mask.
    n_compcor : int
        Number of PCA components per tissue class.

    Returns
    -------
    np.ndarray
        CompCor regression parameters of shape (T, 2 * n_compcor).

    """
    wm_masker = NiftiMasker(mask_img=wm_mask, standardize=True)
    tissue_ts = [wm_masker.fit_transform(fmri_img)]
    if csf_mask is not None:
        csf_masker = NiftiMasker(mask_img=csf_mask, standardize=True)
        tissue_ts.append(csf_masker.fit_transform(fmri_img))

    confounds = []
    for ts in tissue_ts:
        pca = PCA(n_components=min(n_compcor, ts.shape[1]))
        pcs = pca.fit_transform(ts)
        confounds.append(pcs)

    return np.hstack(confounds)


def global_signal(fmri_img, brain_mask):
    """Generate regressor from mean brain signal.

    Parameters
    ----------
    fmri_img : nb.nifti1.Nifti1Image
        Nibabel object containing an fMRI time series.
    brain_mask : nb.nifti1.Nifti1Image
        Nibabel object containing a corresponding brain mask.

    Returns
    -------
    np.ndarray
        Global signal array of shape (T,).

    """
    brain_masker = NiftiMasker(mask_img=brain_mask, standardize=False)
    X = brain_masker.fit_transform(fmri_img)
    return X.mean(axis=1, keepdims=True)


def fd_power(motion, rot_to_mm=50.0, rotations_in_degrees=False):
    """Compute framewise displacement (FD). For FD calculation, rotation (last three
    columns of the motion array) are converted to milimeters by multiplying with head
    radius.

    A scrub mask can be determined by thresholding the FD array, e.g. with a value of
    0.5 mm.

    Parameters
    ----------
    motion : np.ndarray
        Motion parameters of shape (T, 6).
    rot_to_mm : float
        Radius in mm for converting ratations in radians to mm displacement.
    rotations_in_degrees : bool
        If Trye convert rotations to radians first.

    Returns
    -------
    fd : np.ndarray
        FD array of shape (T,).

    """
    motion = motion.copy()
    # Convert rx,ry,rz columns from degrees to radians.
    if rotations_in_degrees:
        motion[:, 3:] = np.deg2rad(motion[:, 3:])

    # Convert rotations to mm (approx)
    motion_mm = motion.copy()
    motion_mm[:, 3:] = motion_mm[:, 3:] * rot_to_mm

    # compute framewise diff
    diff = np.vstack([np.zeros((1, 6)), np.abs(np.diff(motion_mm, axis=0))])
    fd = diff.sum(axis=1)
    return fd


def friston_24(motion):
    """This function generates regressors from motion correction parameters (units do
    not matter here since we want to use these regressors to estimate residuals). 24
    regressors are estimationt: 6 motion parameters + their first temporal derivatives
    + the squares of the motion parameters + the squares of the derivatives.

    Parameters
    ----------
    motion : np.ndarray
        Motion parameters of shape (T, 6).

    Returns
    -------
    confounds : np.ndarray
        Regression parameters of shape (T, 24).

    """
    motion = np.asarray(motion)
    if motion.ndim != 2 or motion.shape[1] != 6:
        raise ValueError("motion must be shape (T,6)")

    # original parameters
    m = motion.copy()

    # temporal derivatives (backwards difference; first row zeros)
    deriv = np.vstack([np.zeros((1, 6)), np.diff(m, axis=0)])

    # squares of the originals and squares of derivatives
    m_sq = m**2
    deriv_sq = deriv**2

    # concatenate in the Friston order: m, deriv, m^2, deriv^2
    confounds = np.hstack([m, deriv, m_sq, deriv_sq])
    return confounds


def get_wm_mask(file_aseg):
    """Get white matter mask from freesurfer aseg file."""
    aseg = nb.load(file_aseg)
    arr_aseg = aseg.get_fdata()
    arr_wm = np.isin(arr_aseg, [2, 41]).astype(np.uint8)
    wm_img = nb.Nifti1Image(arr_wm, affine=aseg.affine)
    return wm_img


def get_csf_mask(file_aseg):
    """Get CSF mask (ventricles) from freesurfer aseg file."""
    aseg = nb.load(file_aseg)
    arr_aseg = aseg.get_fdata()
    arr_csf = np.isin(arr_aseg, [4, 14, 15, 43]).astype(np.uint8)
    csf_img = nb.Nifti1Image(arr_csf, affine=aseg.affine)
    return csf_img
