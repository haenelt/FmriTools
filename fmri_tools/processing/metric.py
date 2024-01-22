# -*- coding: utf-8 -*-
"""Compute metrics from fMRI time series."""

import datetime
import os

import nibabel as nb
import numpy as np
from scipy.stats import zscore

from ..io.filename import get_filename
from ..matlab import MatlabCommand
from ..processing.behavior import get_onset_vols

__all__ = ["compute_psc", "compute_cnr"]


# parameters
USE_LOWPASS = False
CUTOFF_LOWPASS = 0
ORDER_LOWPASS = 0


def compute_psc(
    files_in,
    files_cond,
    TR,
    condition0,
    condition1,
    condition2,
    percent_threshold=50,
    name_output=None,
    name_sess=None,
    use_highpass=False,
    cutoff_highpass=270,
    use_zscore=False,
    skip_vol=2,
    files_outlier=None,
):
    """This scripts calculates the percent signal change (psc) between two conditions of
    a block design for a session consisting of several runs. From the condition file
    which has to be in the SPM compatible *.mat format, time points for both
    experimental blocks and baseline blocks are defined. The psc is computed as the
    difference of the mean between both conditions relative to the the baseline
    condition. The mean psc of the whole session is taken as the average across single
    runs. PSC is computed for both contrast directions. If the outlier input array is
    not empty, outlier volumes are discarded from the analysis. Optionally, the time
    series can be filtered by a lowpass and a highpass filter. The input images should
    be in nifti format.

    Parameters
    ----------
    files_in : list
        List of file names of input time series.
    files_cond : list
        List of file names of condition files in *.mat format.
    TR : float
        Repetition time in s.
    condition0 : str
        Name of baseline condition matching name in condition file.
    condition1 : str
        Name of first experimental condition matching name in condition file.
    condition2 : str
        Name of second experimental condition matching name in condition file.
    percent_threshold : float, optional
        Remove unrealistic high signal values, by default 50.
    name_output : str, optional
        Output name, by default None.
    name_sess : str, optional
        Session name for output file, by default None.
    use_highpass : bool, optional
        Apply highpass filtering, by default False.
    cutoff_highpass : float optional
        Cutoff period for baseline correction in s, by default 270.
    use_zscore : bool, optional
        z-score psc metric, by default False.
    skip_vol : int, optional
        Ignore number of time points at the start of each experimental block, by
        default 2.
    files_outlier : list, optional
        List of file names of regressors of no interest for outlier censoring, by
        default None.
    """
    # get path from first entry
    path_file, _, _ = get_filename(files_in[0])

    # make output folder
    path_output = os.path.join(
        os.path.dirname(os.path.dirname(path_file)), "results", "raw", "native"
    )
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # get image header information
    data_img = nb.load(files_in[0])
    data_img.header["dim"][0] = 3
    data_img.header["dim"][4] = 1
    header = data_img.header
    affine = data_img.affine

    # get image dimension
    dim = data_img.header["dim"][1:4]

    # get outlier dummy array if not outlier input
    if files_outlier is None:
        files_outlier = np.zeros(len(files_in))

    mean_percent_signal1 = np.zeros(dim)
    mean_percent_signal2 = np.zeros(dim)
    for file_, cond_, out_ in zip(files_in, files_cond, files_outlier):
        # get filename
        path_file, name_file, ext_file = get_filename(file_)

        # get condition specific onsets
        onsets0 = get_onset_vols(cond_, out_, condition0, TR, skip_vol)
        onsets1 = get_onset_vols(cond_, out_, condition1, TR, skip_vol)
        onsets2 = get_onset_vols(cond_, out_, condition2, TR, skip_vol)

        # lowpass filter time series
        if USE_LOWPASS:
            matlab = MatlabCommand(
                "ft_lpfilter",
                os.path.join(path_file, f"{name_file}{ext_file}"),
                TR,
                CUTOFF_LOWPASS,
                ORDER_LOWPASS,
            )
            matlab.run()

            # change input to lowpass filtered time series
            name_file = "l" + name_file

        # highpass filter time series
        if use_highpass:
            matlab = MatlabCommand(
                "ft_baseline_correction",
                os.path.join(path_file, f"{name_file}{ext_file}"),
                TR,
                cutoff_highpass,
            )
            matlab.run()

            # change input to highpass filtered time series
            name_file = "b" + name_file

        # open baseline corrected data
        data_img = nb.load(os.path.join(path_file, name_file + ext_file))
        data_array = data_img.get_fdata()

        # sort volumes to conditions
        data_condition0 = data_array[:, :, :, onsets0]
        data_condition1 = data_array[:, :, :, onsets1]
        data_condition2 = data_array[:, :, :, onsets2]

        # z-score
        if use_zscore:
            data_condition0 = zscore(data_condition0, axis=3)
            data_condition1 = zscore(data_condition1, axis=3)
            data_condition2 = zscore(data_condition2, axis=3)

        # mean
        data_condition0_mean = np.mean(data_condition0, axis=3)
        data_condition1_mean = np.mean(data_condition1, axis=3)
        data_condition2_mean = np.mean(data_condition2, axis=3)
        data_condition0_mean[data_condition0_mean == 0] = np.nan

        # percent signal change
        percent_signal1 = (
            (data_condition1_mean - data_condition2_mean) / data_condition0_mean * 100
        )
        percent_signal2 = (
            (data_condition2_mean - data_condition1_mean) / data_condition0_mean * 100
        )

        percent_signal1[np.isnan(percent_signal1)] = 0
        percent_signal2[np.isnan(percent_signal2)] = 0

        # sum volumes for each run
        mean_percent_signal1 += percent_signal1
        mean_percent_signal2 += percent_signal2

    # divide by number of runs
    mean_percent_signal1 /= len(files_in)
    mean_percent_signal2 /= len(files_in)

    # threshold
    if percent_threshold != 0:
        mean_percent_signal1[
            mean_percent_signal1 > percent_threshold
        ] = percent_threshold
        mean_percent_signal1[
            mean_percent_signal1 < -percent_threshold
        ] = -percent_threshold
        mean_percent_signal2[
            mean_percent_signal2 > percent_threshold
        ] = percent_threshold
        mean_percent_signal2[
            mean_percent_signal2 < -percent_threshold
        ] = -percent_threshold

    # name of output files
    if (name_output is not None) and (name_sess is not None):
        fileOUT1 = os.path.join(
            path_output, f"psc_{name_output}_{condition1}_{condition2}_{name_sess}.nii"
        )
        fileOUT2 = os.path.join(
            path_output, f"psc_{name_output}_{condition2}_{condition1}_{name_sess}.nii"
        )
    elif (name_output is not None) and (name_sess is None):
        fileOUT1 = os.path.join(
            path_output,
            f"psc_{name_output}_{condition1}_{condition2}.nii",
        )
        fileOUT2 = os.path.join(
            path_output,
            f"psc_{name_output}_{condition2}_{condition1}.nii",
        )
    elif (name_output is None) and (name_sess is not None):
        fileOUT1 = os.path.join(
            path_output,
            f"psc_{condition1}_{condition2}_{name_sess}.nii",
        )
        fileOUT2 = os.path.join(
            path_output,
            f"psc_{condition2}_{condition1}_{name_sess}.nii",
        )
    else:
        fileOUT1 = os.path.join(path_output, f"psc_{condition1}_{condition2}.nii")
        fileOUT2 = os.path.join(path_output, f"psc_{condition2}_{condition1}.nii")

    # write output
    output = nb.Nifti1Image(mean_percent_signal1, affine, header)
    nb.save(output, fileOUT1)

    output = nb.Nifti1Image(mean_percent_signal2, affine, header)
    nb.save(output, fileOUT2)

    # write log
    file_id = open(os.path.join(path_output, "psc_info.txt"), "a", encoding="utf-8")
    file_id.write(
        f"script executed: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
    )
    file_id.write(f"session: {name_sess}\n")
    file_id.write(f"basename: {name_output}\n")
    file_id.write(f"condition0: {condition0}\n")
    file_id.write(f"condition1: {condition1}\n")
    file_id.write(f"condition2: {condition2}\n")
    if percent_threshold is not None:
        file_id.write(f"percent threshold: {percent_threshold}\n")
    file_id.write(f"TR: {TR}\n")
    file_id.write(f"skip_vol: {skip_vol}\n")
    file_id.write(f"z-score: {use_zscore}\n")
    file_id.write(f"highpass: {use_highpass}\n")
    if use_highpass:
        file_id.write(f"cutoff highpass: {cutoff_highpass}\n")
    if USE_LOWPASS:
        file_id.write(f"lowpass: {USE_LOWPASS}\n")
        file_id.write(f"cutoff lowpass: {CUTOFF_LOWPASS}\n")
        file_id.write(f"order lowpass: {ORDER_LOWPASS}\n")
    file_id.close()


def compute_cnr(
    files_in,
    files_cond,
    TR,
    condition0,
    condition1,
    cnr_threshold,
    name_output=None,
    name_sess=None,
    use_highpass=False,
    cutoff_highpass=270,
    skip_vol=2,
    files_outlier=None,
):
    """This scripts calculates the contrast-to-noise ratio (CNR) in percent from
    functional time series containing task-based activation following a block design.
    The input can be a list of several runs. From the condition file which has to be in
    the SPM compatible *.mat format, time points for both an experimental and a baseline
    condition are extracted. CNR is computed as absolute difference between both
    conditions divided by the standard deviation of the baseline condition. The CNR of
    the whole session is taken as the average across single runs. Similar computations
    of CNR can be found in Scheffler et al. (2016). If the outlier input array is not
    empty, outlier volumes are discarded from the analysis. Optionally, the time series
    can be filtered by a highpass filter. The input images should be in nifti format.
    The script needs an installation of afni.

    Parameters
    ----------
    files_in : list
        List of file names of input time series.
    files_cond : list
        List of file names of condition files in *.mat format.
    TR : float
        Repetition time in s.
    condition0 : str
        Name of baseline condition matching name in condition file.
    condition1 : str
        Name of experimental condition matching name in condition file.
        _description_
    cnr_threshold : float
        Remove unrealistic high signal values.
    name_output : str, optional
        Output name, by default None.
    name_sess : str, optional
        Session name for output file, by default None.
    use_highpass : bool, optional
        Apply highpass filtering, by default False.
    cutoff_highpass : float optional
        Cutoff period for baseline correction in s, by default 270.
    skip_vol : int, optional
        Ignore number of time points at the start of each experimental block, by
        default 2.
    files_outlier : list, optional
        List of file names of regressors of no interest for outlier censoring, by
        default None.
    """
    # get path from first entry
    path_file, _, _ = get_filename(files_in[0])

    # make output folder
    path_output = os.path.join(
        os.path.dirname(os.path.dirname(path_file)), "results", "cnr", "native"
    )
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # get image header information
    data_img = nb.load(files_in[0])
    data_img.header["dim"][0] = 3
    data_img.header["dim"][4] = 1
    header = data_img.header
    affine = data_img.affine

    # get image dimension
    dim = data_img.header["dim"][1:4]

    # get outlier dummy array if not outlier input
    if files_outlier is None:
        files_outlier = np.zeros(len(files_in))

    mean_cnr = np.zeros(dim)
    for file_, cond_, out_ in zip(files_in, files_cond, files_outlier):
        # get filename
        path_file, name_file, ext_file = get_filename(file_)

        # get condition specific onsets
        onsets0 = get_onset_vols(cond_, out_, condition0, TR, skip_vol)
        onsets1 = get_onset_vols(cond_, out_, condition1, TR, skip_vol)

        # highpass filter time series
        if use_highpass:
            matlab = MatlabCommand(
                "ft_baseline_correction",
                os.path.join(path_file, f"{name_file}{ext_file}"),
                TR,
                cutoff_highpass,
            )
            matlab.run()

            # change input to highpass filtered time series
            name_file = "b" + name_file

        # open baseline corrected data
        data_img = nb.load(os.path.join(path_file, f"{name_file}{ext_file}"))
        data_array = data_img.get_fdata()

        # sort volumes to conditions
        data_condition0 = data_array[:, :, :, onsets0]
        data_condition1 = data_array[:, :, :, onsets1]

        # mean
        data_condition0_mean = np.mean(data_condition0, axis=3)
        data_condition1_mean = np.mean(data_condition1, axis=3)
        data_condition0_std = np.std(data_condition0, axis=3)
        data_condition0_std[data_condition0_std == 0] = np.nan

        # percent signal change
        cnr = (
            (np.abs(data_condition1_mean - data_condition0_mean))
            / data_condition0_std
            * 100
        )
        cnr[np.isnan(cnr)] = 0

        # sum volumes for each run
        mean_cnr += cnr

    # divide by number of runs
    mean_cnr /= len(files_in)

    # threshold tsnr
    if cnr_threshold is not None:
        mean_cnr[mean_cnr > cnr_threshold] = cnr_threshold

    # name of output files
    if (name_output is not None) and (name_sess is not None):
        fileOUT = os.path.join(
            path_output, f"cnr_{name_output}_{condition1}_{condition0}_{name_sess}.nii"
        )
    elif (name_output is not None) and (name_sess is None):
        fileOUT = os.path.join(
            path_output,
            f"cnr_{name_output}_{condition1}_{condition0}.nii",
        )
    elif (name_output is None) and (name_sess is not None):
        fileOUT = os.path.join(
            path_output,
            f"cnr_{condition1}_{condition0}_{name_sess}.nii",
        )
    else:
        fileOUT = os.path.join(path_output, f"cnr_{condition1}_{condition0}.nii")

    # write output
    output = nb.Nifti1Image(mean_cnr, affine, header)
    nb.save(output, fileOUT)

    # write log
    file_id = open(os.path.join(path_output, "cnr_info.txt"), "a", encoding="utf-8")
    file_id.write(
        f"script executed: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
    )
    if name_sess:
        file_id.write(f"session: {name_sess}\n")
    if name_output:
        file_id.write(f"basename: {name_output}\n")
    file_id.write(f"condition0: {condition0}\n")
    file_id.write(f"condition1: {condition1}\n")
    if cnr_threshold:
        file_id.write(f"cnr threshold: {cnr_threshold}\n")
    file_id.write(f"TR: {TR}\n")
    file_id.write(f"skip_vol: {skip_vol}\n")
    file_id.write(f"highpass: {use_highpass}\n")
    if use_highpass:
        file_id.write(f"cutoff highpass: {cutoff_highpass}\n")
    file_id.close()
