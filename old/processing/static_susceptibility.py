# -*- coding: utf-8 -*-

import os

import nibabel as nb
import numpy as np

from ..io.filename import get_filename
from ..matlab import MatlabCommand
from ..utils.bias import remove_bias_ants


def static_susceptibility(
    file_in, tr, cutoff_highpass=270, mode="mean", apply_bias=True
):
    """Computation of the temporal mean of an EPI time series. A temporal
    baseline correction and a bias field correction on the resulting temporal
    mean is performed. The resulting image can be used to show locations with
    pronounced static susceptibility effects (veins) in functional time series.

    Parameters
    ----------
    file_in : str
        File name of input time series.
    tr : float
        Volume repetition time (TR) in s.
    cutoff_highpass : float, optional
        Cutoff time in s for baseline correction. Not applied if cutoff is set
        to zero.
    mode : str, optional
        Mode for temporal average (mean or median).
    apply_bias : bool, optional
        Apply bias field correction on temporal mean.

    Returns
    -------
    None.

    """

    if mode not in ["mean", "median"]:
        raise ValueError("Average mode must be either mean or median!")

    # read time series
    data = nb.load(file_in)
    arr = data.get_fdata()

    # output folder structure
    path_file, name_file, ext_file = get_filename(file_in)
    path_output = os.path.join(path_file, "ssw", "native")
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # baseline correction
    if cutoff_highpass != 0:
        matlab = MatlabCommand("ft_baseline_correction", file_in, tr, cutoff_highpass)
        matlab.run()

        # read corrected time series
        arr = nb.load(
            os.path.join(path_file, "b{0}{1}".format(name_file, ext_file))
        ).get_fdata()

    # temporal average
    arr_mean = np.mean(arr, axis=3) if mode == "mean" else np.median(arr, axis=3)

    # write average
    data.header["dim"][0] = 3
    data.header["dim"][4] = 1
    output = nb.Nifti1Image(arr_mean, data.affine, data.header)
    file_out = os.path.join(path_output, "ssw.nii")
    nb.save(output, file_out)

    # bias field correction
    if apply_bias:
        remove_bias_ants(file_out, os.path.join(path_output, "nssw.nii"))
