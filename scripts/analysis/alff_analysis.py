# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import matplotlib.pyplot as plt

# local inputs
#from fmri_tools.analysis.analyze_alff_between_stripes import analyze_alff_between_stripes
from fmri_tools.analysis.analyze_alff_between_conditions import analyze_alff_between_conditions


"""
Comparison of resting-state data within and between V2 stripes.

created by Daniel Haenelt
Date created: 11-03-2019
Last modified: 12-10-2020
"""

# input paths
input_label = '/media/haenelt/6CFE-5CB9/p6_resting_state/v2/lh.v2.label'
input_color = '/data/pt_01880/V2STRIPES/p6/colour/results/spmT/surf/lh.spmT_colour_bw_GE_EPI1_def_layer5.mgh'
input_disparity = '/data/pt_01880/V2STRIPES/p6/disparity/results/spmT/surf/lh.spmT_disparity_control_GE_EPI1_def_layer5.mgh'
input_rest = '/data/pt_01880/V2STRIPES/p6/resting_state/alff_biopac/surf/lh.falff_z_def-img.nii_layer5_def.mgh'
path_output = '/data/pt_01880'

# parameters
min_contrast = 2.7
nvert = 1000

# do not edit below

# compare alff within and between stripes
#rest1, rest2, t, p, _ = analyze_alff_between_stripes(input_label, input_disparity, input_rest, 
#                                                     min_contrast, nvert)

rest1, rest2, t, p, _ = analyze_alff_between_conditions(input_label, 
                                                        input_color, 
                                                        input_disparity, 
                                                        input_rest, 
                                                        min_contrast, 
                                                        nvert)

# boxplot
data = [rest1, rest2]
fig = plt.figure()
plt.boxplot(data, showfliers=False)
plt.rcParams["font.size"]= 10
plt.xticks([1,2],["color","disparity"])
plt.ylabel("z-score (fALFF)")
plt.figtext(0.27, 0.8, "t = "+str(np.round(t,3))+", p = "+str(np.round(p, 3)), horizontalalignment="center")
plt.savefig(os.path.join(path_output,"falff_biopac_condition.png"),dpi=100)
plt.show()
