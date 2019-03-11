"""
Comparison of resting-state data within and between V2 stripes.

created by Daniel Haenelt
Date created: 11-03-2019
Last modified: 11-03-2019
"""
import os
import numpy as np
import nibabel as nb
import random
import matplotlib.pyplot as plt
from nibabel.freesurfer.io import read_label
from scipy.stats import ttest_ind, levene

# input paths
input_label = '/data/pt_01880/lh.v2.label'
input_contrast = '/data/pt_01880/V2STRIPES/p6/colour/results/spmT/surf/lh.spmT_colour_bw_GE_EPI1_def_layer5.mgh'
input_rest = '/data/pt_01880/V2STRIPES/p6/resting_state/alff/surf/lh.alff_z_def-img.nii_layer5_def.mgh'
path_output = '/data/pt_01880'

# parameters
min_contrast = 2.7
nvert = 1000

""" do not edit below """

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

# get resting-state data in mask
rest_pos = rest.copy()
rest_neg = rest.copy()

rest_pos = rest_pos * contrast_pos
rest_pos = rest_pos[~np.isnan(rest_pos)]

rest_neg = rest_neg * contrast_neg
rest_neg = rest_neg[~np.isnan(rest_neg)]

# select random number of vertices
label_shuffled_pos = random.sample(range(0,len(rest_pos)), nvert)
label_shuffled_neg = random.sample(range(0,len(rest_neg)), nvert)

rest_pos = rest_pos[label_shuffled_pos]
rest_neg = rest_neg[label_shuffled_neg]

# independent samples t-test
# Levene's test is run to check for equal variances. If variances are not equal, Welch's t-test
# is performed.
_, p_levene = levene(rest_pos, rest_neg)
if p_levene < 0.05:
    t, p = ttest_ind(rest_pos, rest_neg, equal_var=False)
else:
    t, p = ttest_ind(rest_pos, rest_neg, equal_var=True)

# boxplot
data = [rest_pos, rest_neg]
fig = plt.figure()
plt.boxplot(data)
plt.rcParams["font.size"]= 10
plt.xticks([1,2],["Within stripes","Between stripes"])
plt.ylabel("z-score (ALFF)")
plt.figtext(0.27, 0.8, "t = "+str(np.round(t,3))+", p = "+str(np.round(p, 3)), horizontalalignment="center")
plt.savefig(os.path.join(path_output,"alff_stripes.png"),dpi=100)
plt.show()