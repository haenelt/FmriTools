# -*- coding: utf-8 -*-
"""
Analysis of behavioural data

Analysis of the behavioural performance of the dummy fixation task during
scanning. The script essentially performs the same calculation as directly done
automatically after single functional runs.

"""

import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
from scipy.io import loadmat

file_in = [
    "/data/pt_01880/Experiment1_ODC/p1/odc/GE_EPI1/Run_1/logfiles/p1_GE_EPI1_Run1_odc.mat",
]

# do not edit below

# font parameters for plots
rc("font", **{"family": "serif", "serif": ["Palatino"]})
rc("text", usetex=True)

for f in file_in:
    # get output folder
    path_output = os.path.join(os.path.dirname(os.path.dirname(f)), "behavior")
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # get mat-file
    FixationData = loadmat(f)["FixationData"]

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
    fileID = open(os.path.join(path_output, "fixation_task_summary.txt"), "w")
    fileID.write("Number of changes: %i\n" % change)
    fileID.write("Number of responses: %i\n" % response)
    fileID.write("Number of hits: %i\n" % change_hit)
    fileID.write("Number of misses: %i\n" % change_miss)
    fileID.write("Error rate: %.2f %%\n" % error_rate)
    fileID.write("Mean RT: %.2f ms\n" % mean_rt)
    fileID.write("Corresponding SD: %.2f ms" % std_rt)
    fileID.close()

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
