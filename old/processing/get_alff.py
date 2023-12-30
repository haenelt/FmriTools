# -*- coding: utf-8 -*-

import os

import nibabel as nb
from nipype.interfaces.afni.preprocess import Bandpass
from nipype.interfaces.afni.utils import Calc, TStat
from scipy.stats import zscore

from ..io.filename import get_filename


def get_alff(file_in, TR, path_output, hp_freq=0.01, lp_freq=0.08, cleanup=True):
    """Get ALFF.

    This function calculates ALFF and fALFF from a preprocessed (motion
    correction, nuisance regression, etc.) resting-state time series. ALFF is
    computed by bandpass filtering the time series and computing the voxel-wise
    standard deviation of the filtered time series. fALFF is computed by
    dividing ALFF by the voxel-wise standard deviation of the unfiltered time
    series. Additionally, ALFF and fALFF are expressed in z-score. This function
    follows the script found in [1].

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

    Returns
    -------
    None.

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
    bandpass = Bandpass()
    bandpass.inputs.in_file = file_in
    bandpass.inputs.apply_inverse = hp_freq
    bandpass.inputs.apply = lp_freq
    bandpass.inputs.tr = TR
    bandpass.inputs.outputtype = "NIFTI"
    bandpass.inputs.out_file = os.path.join(path_output, file + "_filtered.nii")
    bandpass.run()

    # standard deviation over frequency
    stddev_filtered = TStat()
    stddev_filtered.inputs.in_file = os.path.join(path_output, file + "_filtered.nii")
    stddev_filtered.inputs.args = "-stdev"
    stddev_filtered.inputs.outputtype = "NIFTI"
    stddev_filtered.inputs.out_file = os.path.join(path_output, "alff.nii")
    stddev_filtered.run()

    # standard deviation of the unfiltered nuisance corrected image
    stddev_unfiltered = TStat()
    stddev_unfiltered.inputs.in_file = file_in
    stddev_unfiltered.inputs.args = "-stdev"
    stddev_unfiltered.inputs.outputtype = "NIFTI"
    stddev_unfiltered.inputs.out_file = os.path.join(path_output, "temp.nii")
    stddev_unfiltered.run()

    # falff calculations
    falff = Calc()
    falff.inputs.in_file_a = os.path.join(path_output, "alff.nii")
    falff.inputs.in_file_b = os.path.join(path_output, "temp.nii")
    falff.inputs.args = "-float"
    falff.inputs.expr = "(1.0*a)/(1.0*b)"
    falff.inputs.outputtype = "NIFTI"
    falff.inputs.out_file = os.path.join(path_output, "falff.nii")
    falff.run()

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
