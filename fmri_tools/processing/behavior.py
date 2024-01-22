# -*- coding: utf-8 -*-
"""Behavioral measures."""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
from scipy.io import loadmat

__all__ = ["analysis_fixation", "get_onset_vols"]

# font parameters for plots
rc("font", **{"family": "serif", "serif": ["Palatino"]})
rc("text", usetex=True)


# mat condition file
def analysis_fixation(file_in):
    """Analysis of the behavioural performance of the dummy fixation task during
    scanning. The script essentially performs the same calculation as directly done
    automatically after single functional runs.

    Parameters
    ----------
    file_in : str
        File name of condition file in mat format containing the key FixationData.
    """
    # get output folder
    path_output = os.path.join(os.path.dirname(os.path.dirname(file_in)), "behavior")
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # get mat-file
    FixationData = loadmat(file_in)["FixationData"]

    # number of responses
    response = 0
    change = 0
    for j in range(len(FixationData)):
        if FixationData[j, 0] == 3:
            response += 1
        else:
            change += 1

    # number of hits, misses and rt
    change_miss = 0
    change_hit = 0
    rt = []
    for j in range(len(FixationData) - 1):
        if FixationData[j, 0] != 3 and FixationData[j + 1, 0] != 3:
            change_miss += 1
        elif FixationData[j, 0] != 3 and FixationData[j + 1, 0] == 3:
            change_hit += 1
            rt.append((FixationData[j + 1, 1] - FixationData[j, 1]) * 1000)

    # add a miss if the run does not end with a response
    if FixationData[-1, 0] != 3:
        change_miss += 1

    # compute error rate
    error_rate = change_miss / change * 100

    # reaction time mean and std
    mean_rt = np.mean(rt)
    std_rt = np.std(rt)

    # output
    file_id = open(
        os.path.join(path_output, "fixation_task_summary.txt"), "w", encoding="utf-8"
    )
    file_id.write(f"Number of changes: {change}\n")
    file_id.write(f"Number of responses: {response}\n")
    file_id.write(f"Number of hits: {change_hit}\n")
    file_id.write(f"Number of misses: {change_miss}\n")
    file_id.write(f"Error rate: {error_rate:.2f} %\n")
    file_id.write(f"Mean RT: {mean_rt:.2f} ms\n")
    file_id.write(f"Corresponding SD: {std_rt:.2f} ms")
    file_id.close()

    # hist plot
    fig, ax = plt.subplots()
    ax.hist(rt)
    ax.set_xlabel("RT in ms")
    ax.set_ylabel("Number of responses")
    ax.set_title("Dummy fixation task reaction times")
    fig.savefig(
        os.path.join(path_output, "fixation_task_hist.png"),
        format="png",
        bbox_inches="tight",
    )


def get_onset_vols(cond_input, outlier_input, name_condition, TR, skip_vol):
    """This function returns all the volume indices corresponding to an experimental
    condition in from a block design.

    Parameters
    ----------
    cond_input : str
        Block design condition *.mat file.
    outlier_input : str
        Regressor of no interest *.txt file.
    name_condition : str
        Name of first condition.
    TR : float
        Repetition time in s.
    skip_vol : int
        Number of skipped time point of each block.

    Returns
    -------
    onsets : ndarray
        sorted volumes of experimental condition.

    """
    # load condition file
    cond = loadmat(cond_input)

    # get condition information
    names = np.concatenate(np.concatenate(cond["names"]))
    onsets = np.concatenate(cond["onsets"])
    durations = np.concatenate(np.concatenate(np.concatenate(cond["durations"])))

    # check if condition names exist
    if name_condition not in names:
        sys.exit("The condition is not found in the condition_file")

    # index of condition
    ind = np.where(name_condition == names)[0][0]

    # onsets for both conditions
    onsets = np.round(onsets[ind][0] / TR + skip_vol).astype(int)

    # durations for both conditions
    durations = np.round(durations[ind] / TR - skip_vol).astype(int)

    # sort all volumes to be considered in both conditions
    temp = onsets.copy()
    for i in range(durations - 1):
        onsets = np.append(onsets, temp + i + 1)
    onsets = np.sort(onsets)

    # remove outlier volumes
    if outlier_input:
        # load outlier regressor
        outlier_regressor = np.loadtxt(outlier_input)
        outlier_regressor = np.where(outlier_regressor == 1)[0]

        # look for outliers in onset arrays
        for i, _ in enumerate(onsets):
            if np.any(onsets[i] == outlier_regressor):
                onsets[i] = -1

        # remove outliers
        onsets = onsets[onsets != -1]

    return onsets
