# -*- coding: utf-8 -*-
"""Hemodynamic response function."""

import numpy as np
from scipy.stats import gamma

__all__ = ["hrf_spm"]


def hrf_spm(TR, p=(6, 16, 1, 1, 6, 0, 32)):
    """An implementation of spm_hrf.m from the SPM distribution. The code is taken from
    [1]_.

    Parameters
    ----------
    TR : float
        Repetition time in s.
    p : tuple, optional
        Parameters of the two gamma functions.

    Default settings of p in seconds:
    p[0] - delay of response (relative to onset)         6
    p[1] - delay of undershoot (relative to onset)      16
    p[2] - dispersion of response                        1
    p[3] - dispersion of undershoot                      1
    p[4] - ratio of response to undershoot               6
    p[5] - onset (seconds)                               0
    p[6] - length of kernel (seconds)                   32

    Returns
    -------
    hrf : ndarray
        The HRF.

    References
    -------
    .. [1] https://github.com/poldracklab/poldracklab-base/blob/master/fmri/spm_hrf.py

    """
    p = [float(x) for x in p]
    fMRI_T = 16.0
    TR = float(TR)
    dt = TR / fMRI_T

    u = np.arange(p[6] / dt + 1) - p[5] / dt
    hrf = (
        gamma.pdf(u, p[0] / p[2], scale=1.0 / (dt / p[2]))
        - gamma.pdf(u, p[1] / p[3], scale=1.0 / (dt / p[3])) / p[4]
    )

    good_pts = np.array(range(int(p[6] / TR))) * fMRI_T
    good_pts = good_pts.astype(int)  # das habe ich eingefuegt
    hrf = hrf[list(good_pts)]
    # hrf = hrf([0:(p(7)/RT)]*fMRI_T + 1);
    hrf = hrf / np.sum(hrf)

    return hrf
