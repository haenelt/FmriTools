# -*- coding: utf-8 -*-
"""Estimate scan-to-scan repeatability."""

import os
import random

import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
from nibabel.freesurfer.io import read_label, read_morph_data
from scipy.stats import pearsonr


def estimate_reliability(
    file_sess1,
    file_sess2,
    file_label,
    path_output,
    basename_output,
    frac=0.25,
    niter=1000,
):
    """Reliability measure between sessions. Data within a defined label is compared
    vertex-wise using Pearson correlation. Only a fraction of vertices is used to
    account for spatial covariance of nearby vertices. The p-value is estimated using a
    permutation analysis. A scatter and a Bland-Altman plot is created.

    Parameters
    ----------
    file_sess1 : str
        MGH overlay of first session.
    file_sess2 : str
        MGH overlay of second session.
    file_label : str
        Freesurfer label file.
    path_output : str
        Path where output is written.
    basename_output : _type_
        Basename of output files.
    frac : float, optional
        Fraction of considered vertices in label, by default 0.25.
    niter : int, optional
        Number of iterations, by default 1000.
    """
    # correlation plot labels
    corr_title = ""
    corr_x_label = "session 1"
    corr_y_label = "session 2"

    # load data
    label = read_label(file_label).tolist()

    # if input file extension is not *.mgh, interprete as morphological file
    if os.path.splitext(os.path.basename(file_sess1))[1] == ".mgh":
        sess1 = np.squeeze(nb.load(file_sess1).get_fdata())
    else:
        sess1 = np.squeeze(read_morph_data(file_sess1))

    # if input file extension is not *.mgh, interprete as morphological file
    if os.path.splitext(os.path.basename(file_sess2))[1] == ".mgh":
        sess2 = np.squeeze(nb.load(file_sess2).get_fdata())
    else:
        sess2 = np.squeeze(read_morph_data(file_sess2))

    # get the amount of data points
    ndata = np.round(frac * len(label)).astype(int)

    # randomly select ndata points in sess1 and sess2
    label_shuffled = random.sample(label, len(label))
    label_shuffled = label_shuffled[0:ndata]

    # get correlation between sessions
    x = sess1[label_shuffled]
    y = sess2[label_shuffled]
    r = pearsonr(x, y)

    # Bland Altman plot
    mean = np.mean([x, y], axis=0)
    diff = x - y  # difference between data1 and data2
    md = np.mean(diff)  # mean of the difference
    sd = np.std(diff, axis=0)  # Standard deviation of the difference

    plt.scatter(mean, diff)
    plt.axhline(md, color="gray", linestyle="--")
    plt.axhline(md + 1.96 * sd, color="gray", linestyle="--")
    plt.axhline(md - 1.96 * sd, color="gray", linestyle="--")
    plt.title("Bland-Altman Plot: " + os.path.basename(file_label))
    plt.xlabel("Average of 2 sessions")
    plt.ylabel("Difference between 2 sessions")
    plt.savefig(os.path.join(path_output, basename_output + "_bland_altman.png"))

    # compare correlation coefficient to change level (permutation)
    null_dist = []
    for i in range(niter):
        # permute data points in sess2
        label_shuffled = random.sample(label, len(label))
        label_shuffled = label_shuffled[0:ndata]

        y_shuffle = sess2[label_shuffled]
        null_dist.append(pearsonr(x, y_shuffle)[0])

    p = np.array(null_dist) > r[0]
    p = p[p is True]
    p = len(p) / niter

    # print results
    print("ROI: " + os.path.basename(file_label))
    print("Number of points: " + str(ndata))
    print("Number of iterations: " + str(niter))
    print("Correlation coefficient: " + str(r[0]))
    print("p-value: " + str(p))

    # scatter plot with regression line
    plt.figure()
    plt.scatter(x, y)
    plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)), "r")
    plt.title(corr_title)
    plt.xlabel(corr_x_label)
    plt.ylabel(corr_y_label)
    plt.figtext(
        0.27,
        0.8,
        "r = " + str(np.round(r[0], 3)) + ", p = " + str(np.round(p, 3)),
        horizontalalignment="center",
    )
    plt.savefig(os.path.join(path_output, basename_output + ".png"))
    plt.show()
