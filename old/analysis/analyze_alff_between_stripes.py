# -*- coding: utf-8 -*-

import random

import nibabel as nb
import numpy as np
from nibabel.freesurfer.io import read_label
from scipy.stats import levene, ttest_ind


def analyze_alff_between_stripes(
    input_label, input_contrast, input_rest, min_contrast, nvert
):
    """Analyze ALFF between stripes.

    Comparison of resting-state data within and between V2 stripes. Mask stripes
    within a ROI from a label with a thresholded contrast and select randomly
    <nvert> vertices either within or between stripes. An independent samples
    t-test is computed (or Welch's test if Levene's test is significant).

    Parameters
    ----------
    input_label : str
         Input label file to define the region of interest.
    input_contrast : str
        Input contrast data for masking.
    input_rest : str
        Input resting-state data (alff or falff).
    min_contrast : float
        Minimum contrast for masking (t-score).
    nvert : int
        Number of selected vertices.

    Returns
    -------
    rest_pos : ndarray
        Resting-state data within stripes.
    rest_neg : ndarray
        Resting-state data between stripes.
    t : float
        t-score from independent samples t-test.
    p : float
        p-value from independent samples t-test.
    p_levene : float
        p-value from Levene's test.

    """

    # load data
    label = read_label(input_label)
    contrast = np.squeeze(nb.load(input_contrast).get_fdata())
    rest = np.squeeze(nb.load(input_rest).get_fdata())

    # mask data
    contrast_pos = contrast.copy()
    contrast_neg = contrast.copy()

    contrast_pos[contrast < min_contrast] = np.NaN
    contrast_neg[contrast > -min_contrast] = np.NaN

    contrast_pos[~label] = np.NaN
    contrast_neg[~label] = np.NaN
    contrast_pos[~np.isnan(contrast_pos)] = 1
    contrast_neg[~np.isnan(contrast_neg)] = 1

    # get resting-state data in mask
    rest_pos = rest.copy()
    rest_neg = rest.copy()

    rest_pos = rest_pos * contrast_pos
    rest_pos = rest_pos[~np.isnan(rest_pos)]
    rest_pos = rest_pos[rest_pos != np.min(rest_pos)]

    rest_neg = rest_neg * contrast_neg
    rest_neg = rest_neg[~np.isnan(rest_neg)]
    rest_neg = rest_neg[rest_neg != np.min(rest_neg)]

    # select random number of vertices
    label_shuffled_pos = random.sample(range(0, len(rest_pos)), nvert)
    label_shuffled_neg = random.sample(range(0, len(rest_neg)), nvert)

    rest_pos = rest_pos[label_shuffled_pos]
    rest_neg = rest_neg[label_shuffled_neg]

    # independent samples t-test
    # Levene's test is run to check for equal variances. If variances are not
    # equal, Welch's t-test is performed.
    _, p_levene = levene(rest_pos, rest_neg)
    if p_levene < 0.05:
        t, p = ttest_ind(rest_pos, rest_neg, equal_var=False)
    else:
        t, p = ttest_ind(rest_pos, rest_neg, equal_var=True)

    return rest_pos, rest_neg, t, p, p_levene
