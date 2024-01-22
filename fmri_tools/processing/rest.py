# -*- coding: utf-8 -*-
"""Resting-state utilities."""

import os
import subprocess

import nibabel as nb
from scipy.stats import zscore

from ..io.filename import get_filename
from ..preprocessing.timeseries import bandpass_afni

__all__ = ["get_alff"]


def get_alff(file_in, TR, path_output, hp_freq=0.01, lp_freq=0.08, cleanup=True):
    """This function calculates ALFF and fALFF from a preprocessed (motion correction,
    nuisance regression, etc.) resting-state time series. ALFF is computed by bandpass
    filtering the time series and computing the voxel-wise standard deviation of the
    filtered time series. fALFF is computed by dividing ALFF by the voxel-wise standard
    deviation of the unfiltered time series. Additionally, ALFF and fALFF are expressed
    in z-score. This function follows the script found in [1].

    Parameters
    ----------
    file_in : str
        Input time series.
    TR : float
        Repetition time in s.
    path_output : str
        Path where output is written.
    hp_freq : float, optional
        Highpass cutoff frequency in Hz. The default is 0.01.
    lp_freq : float, optional
        Lowpass cutoff frequency in Hz. The default is 0.08.
    cleanup : bool, optional
        Delete intermediate files. The default is True.

    References
    -------
    .. [1] https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/alff/alff.py

    """
    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # get path and filename
    _, file, _ = get_filename(file_in)

    # filtering
    bandpass_afni(
        file_in, os.path.join(path_output, file + "_filtered.nii"), TR, lp_freq, hp_freq
    )

    # standard deviation over frequency
    _std_afni(
        os.path.join(path_output, file + "_filtered.nii"),
        os.path.join(path_output, "alff.nii"),
    )

    # standard deviation of the unfiltered nuisance corrected image
    _std_afni(file_in, os.path.join(path_output, "temp.nii"))

    # falff calculations
    _divide_afni(
        os.path.join(path_output, "alff.nii"),
        os.path.join(path_output, "temp.nii"),
        os.path.join(path_output, "falff.nii"),
    )

    # alff in z-score
    alff_img = nb.load(os.path.join(path_output, "alff.nii"))
    alff_array = alff_img.get_fdata()
    alff_array = zscore(alff_array, axis=None)

    output = nb.Nifti1Image(alff_array, alff_img.affine, alff_img.header)
    nb.save(output, os.path.join(path_output, "alff_z.nii"))

    # falff in z-score
    falff_img = nb.load(os.path.join(path_output, "falff.nii"))
    falff_array = falff_img.get_fdata()
    falff_array = zscore(falff_array, axis=None)

    output = nb.Nifti1Image(falff_array, falff_img.affine, falff_img.header)
    nb.save(output, os.path.join(path_output, "falff_z.nii"))

    # cleanup
    if cleanup:
        os.remove(os.path.join(path_output, "temp.nii"))
        os.remove(os.path.join(path_output, file + "_filtered.nii"))


def _divide_afni(file_in1, file_in2, file_out):
    """Divide two images using AFNI.

    Parameters
    ----------
    file_in1 : str
        File name of first input file.
    file_in2 : str
        File name of second input file.
    file_out : str
        File name of output file.
    """
    command = "3dcalc"
    command += f" -a {file_in1}"
    command += f" -b {file_in2}"
    command += " -expr '(1.0*a)/(1.0*b)'"
    command += " -float"
    command += f" -prefix {file_out}"

    print("Execute: " + command)
    try:
        subprocess.run([command], shell=True, check=False)
    except subprocess.CalledProcessError:
        print("Execuation failed!")


def _std_afni(file_in, file_out):
    """Compute voxel-wise standard deviation using AFNI.

    Parameters
    ----------
    file_in : str
        File name of input file.
    file_out : str
        File name of output file.
    """
    command = "3dTstat"
    command += " -stdev"
    command += f" -prefix {file_out}"
    command += f" {file_in}"

    print("Execute: " + command)
    try:
        subprocess.run([command], shell=True, check=False)
    except subprocess.CalledProcessError:
        print("Execuation failed!")
