# -*- coding: utf-8 -*-
"""Tools for smoothing 3D nifti data."""

from pathlib import Path
import numpy as np
import nibabel as nb
from dipy.denoise.nlmeans import nlmeans
from dipy.denoise.noise_estimate import estimate_sigma

__all__ = ["nlm"]


def nlm(file_in, file_out, n=64, scale=1.0):
    """Apply edge-preserving smoothing to nifti data by non-local means filtering.

    Parameters
    ----------
    file_in: str
        File name of input file.
    file_out: str
        File name of output file.
    n: int, optional
        Number of coils of the receiver array.
    scale: float, optional
        Scale factor for noise standard deviation estimation.
    """
    data = nb.load(file_in)
    arr = data.get_fdata()
    sigma = estimate_sigma(arr, N=n)  # estimate noise standard deviation
    sigma *= scale
    arr_smoothed = nlmeans(arr, sigma=sigma)
    # write output with integer datatype
    header = data.header
    header.set_data_dtype(np.float32)
    arr_smoothed = arr_smoothed.astype(np.float32)
    output = nb.Nifti1Image(arr_smoothed, data.affine, header)
    Path(file_out).parent.mkdir(exist_ok=True, parents=True)
    nb.save(output, file_out)
