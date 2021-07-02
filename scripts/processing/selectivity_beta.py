# -*- coding: utf-8 -*-
"""
Selectivity metric a la Olman

Measures for BOLD sensitivity and selectivity between two orthogonal conditions
are computed. The metrics are adapted from Olman et al. 2018. Sensitivity is
defined as the average response (betas) to the visual stimulus (sensitivity =
(|cond1| + |cond2|) / 2). This metric is used to assess responsivity of the
gray matter throughout cortical depth. Selectivity is defined as the average
absolute difference of both conditions (betas) divided by the average response.

"""

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb

# input data
input_beta_cond1 = "/data/pt_01880/Experiment1_ODC/p3/odc/results/con/native/con_left_GE_EPI2_test.nii"
input_beta_cond2 = "/data/pt_01880/Experiment1_ODC/p3/odc/results/con/native/con_right_GE_EPI2_test.nii"

# parameters
condition1 = "left"
condition2 = "right"
name_sess = "GE_EPI2"
name_output = ""
path_output = "/data/pt_01880/"

# do not edit below

# output folder is taken from the first entry of the input list
path_output = os.path.join(path_output, "results", "selectivity", "native")
if not os.path.exists(path_output):
    os.makedirs(path_output)

# load input
cond1 = nb.load(input_beta_cond1)
cond2 = nb.load(input_beta_cond2)

cond1_img = cond1.get_fdata()
cond2_img = cond2.get_fdata()

# get activtation metric
sensitivity = (np.abs(cond1_img) + np.abs(cond2_img)) / 2
selectivity = np.abs(cond1_img - cond2_img) / 2 / sensitivity

# name of output files
if len(name_output) and len(name_sess):
    fileOUT1 = os.path.join(path_output,
                            "sensitivity_beta_" +
                            name_output + "_" +
                            condition1 + "_" +
                            condition2 + "_" +
                            name_sess + ".nii")
    fileOUT2 = os.path.join(path_output,
                            "selectivity_beta_" +
                            name_output + "_" +
                            condition1 + "_" +
                            condition2 + "_" +
                            name_sess + ".nii")
elif len(name_output) and not len(name_sess):
    fileOUT1 = os.path.join(path_output,
                            "sensitivity_beta_" +
                            name_output + "_" +
                            condition1 + "_" +
                            condition2 + ".nii")
    fileOUT2 = os.path.join(path_output,
                            "selectivity_beta_" +
                            name_output + "_" +
                            condition1 + "_" +
                            condition2 + ".nii")
elif not len(name_output) and len(name_sess):
    fileOUT1 = os.path.join(path_output,
                            "sensitivity_beta_" +
                            condition1 + "_" +
                            condition2 + "_" +
                            name_sess + ".nii")
    fileOUT2 = os.path.join(path_output,
                            "selectivity_beta_" +
                            condition1 + "_" +
                            condition2 + "_" +
                            name_sess + ".nii")
else:
    fileOUT1 = os.path.join(path_output,
                            "sensitivity_beta_" +
                            condition1 + "_" +
                            condition2 + ".nii")
    fileOUT2 = os.path.join(path_output,
                            "selectivity_beta_" +
                            condition1 + "_" +
                            condition2 + ".nii")

# write output
output1 = nb.Nifti1Image(sensitivity, cond1.affine, cond1.header)
output2 = nb.Nifti1Image(selectivity, cond1.affine, cond1.header)
nb.save(output1, fileOUT1)
nb.save(output2, fileOUT2)
