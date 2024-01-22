# -*- coding: utf-8 -*-
"""Estimate bias metrics in fMRI data."""

import os

import nibabel as nb
import numpy as np

from ..io.filename import get_filename
from ..matlab import MatlabCommand
from ..segmentation.vol import remove_bias_ants

__all__ = ["estimate_pv", "static_susceptibility"]


def estimate_pv(input_target, input_border, path_output, name_output):
    """This function estimates the partial volume contribution in each image voxel of a
    target image from an upsampled binary image depicting the GM/WM or GM/CSF border.
    Partial voluming is estimated by downsampling the binary image using a
    moving-average like algorithm and calculating the ratio of both binary elements
    within each target voxel.

    Parameters
    ----------
    input_target : str
        Target space for which partial voluming is estimated.
    input_border : str
        Upsampled binary border depicting the high-resolution tissue border.
    path_output : str
        Path where output is saved.
    name_output : str
        Basename of output image.
    """
    # load data
    target = nb.load(input_target)
    border = nb.load(input_border)
    border_array = border.get_fdata()

    # downsampling parameters
    matrix_up = border.header["dim"][1:4]
    matrix_down = target.header["dim"][1:4]
    magn_factor = matrix_up / matrix_down

    # sampling steps
    p = np.arange(0, matrix_down[0] * magn_factor[0], magn_factor[0])
    q = np.arange(0, matrix_down[1] * magn_factor[1], magn_factor[1])
    r = np.arange(0, matrix_down[2] * magn_factor[2], magn_factor[2])

    # initialise matrices
    M = np.zeros(target.header["dim"][1:4])
    start_weight = np.zeros(3)
    end_weight = np.zeros(3)

    # downsampling using a moving-average like algorithm
    for i, _ in enumerate(p):
        # weightings (start_weight)
        if np.mod(p[i], 1) != 0:
            start_weight[0] = 1 - np.mod(p[i], 1)
        else:
            start_weight[0] = 1

        # weightings (end_weight)
        if i == len(p) - 1:
            end_weight[0] = 1
        else:
            if np.mod(p[i + 1], 1) != 0:
                end_weight[0] = np.mod(p[i + 1], 1)
            else:
                end_weight[0] = 1

        for j, _ in enumerate(q):
            # weightings (start_weight)
            if np.mod(q[j], 1) != 0:
                start_weight[1] = 1 - np.mod(q[j], 1)
            else:
                start_weight[1] = 1

            # weightings (end_weight)
            if j == len(q) - 1:
                end_weight[1] = 1
            else:
                if np.mod(q[j + 1], 1) != 0:
                    end_weight[1] = np.mod(q[j + 1], 1)
                else:
                    end_weight[1] = 1

            for k, _ in enumerate(r):
                # weightings (start_weight)
                if np.mod(r[k], 1) != 0:
                    start_weight[2] = 1 - np.mod(r[k], 1)
                else:
                    start_weight[2] = 1

                # weightings (end_weight)
                if k == len(r) - 1:
                    end_weight[2] = 1
                else:
                    if np.mod(r[k + 1], 1) != 0:
                        end_weight[2] = np.mod(r[k + 1], 1)
                    else:
                        end_weight[2] = 1

                # p-coordinates
                if i == len(p) - 1:
                    x = np.arange(np.round(p[-1] - 0.5), matrix_up[0], 1).astype(int)
                else:
                    x = np.arange(
                        np.round(p[i] - 0.5), np.round(p[i + 1] - 0.5), 1
                    ).astype(int)

                # q-coordinates
                if j == len(q) - 1:
                    y = np.arange(np.round(q[-1] - 0.5), matrix_up[1], 1).astype(int)
                else:
                    y = np.arange(
                        np.round(q[j] - 0.5), np.round(q[j + 1] - 0.5), 1
                    ).astype(int)

                # r-coordinates
                if k == len(r) - 1:
                    z = np.arange(np.round(r[-1] - 0.5), matrix_up[2], 1).astype(int)
                else:
                    z = np.arange(
                        np.round(r[k] - 0.5), np.round(r[k + 1] - 0.5), 1
                    ).astype(int)

                # define submatrix
                X = border_array[x, :, :]
                X = X[:, y, :]
                X = X[:, :, z]

                # apply weightings
                X[0, :, :] = start_weight[0] * X[0, :, :]
                X[:, 0, :] = start_weight[1] * X[:, 0, :]
                X[:, :, 0] = start_weight[2] * X[:, :, 0]
                X[-1, :, :] = end_weight[0] * X[-1, :, :]
                X[:, -1, :] = end_weight[1] * X[:, -1, :]
                X[:, :, -1] = end_weight[2] * X[:, :, -1]

                X = np.reshape(X, np.size(X))
                X = X[~np.isnan(X)]

                # voxel value in target space
                M[i, j, k] = np.sum(X) / len(X)

    # save data
    output = nb.Nifti1Image(M, target.affine, target.header)
    nb.save(output, os.path.join(path_output, name_output + "_pve.nii"))


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
