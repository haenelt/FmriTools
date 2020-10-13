# -*- coding: utf-8 -*-

# local inputs
from fmri_tools.preprocessing.slice_timing_correction import slice_timing_correction


"""
Temporal correction

This scripts performs slice timing correction and temporal regridding on a nifti 
time series.

The script needs an installation of afni.

created by Daniel Haenelt
Date created: 11-03-2020
Last modified: 13-10-2020
"""

# input
input = ["/data/pt_01880/Experiment1_ODC/p4/retinotopy/pol_anticlock/data.nii",
         "/data/pt_01880/Experiment1_ODC/p4/retinotopy/pol_clock/data.nii",
         "/data/pt_01880/Experiment1_ODC/p4/retinotopy/ecc_expanding/data.nii",
         "/data/pt_01880/Experiment1_ODC/p4/retinotopy/ecc_contracting/data.nii",
         "/data/pt_01880/Experiment1_ODC/p4/retinotopy2/pol_anticlock/data.nii",
         "/data/pt_01880/Experiment1_ODC/p4/retinotopy2/pol_clock/data.nii",
         "/data/pt_01880/Experiment1_ODC/p4/retinotopy2/ecc_expanding/data.nii",
         "/data/pt_01880/Experiment1_ODC/p4/retinotopy2/ecc_contracting/data.nii",
         ]
TR_old = 2.0
TR_new = 2.0
order = "descending"

# do not edit below

for i in range(len(input)):
    slice_timing_correction(input[i], TR_old, TR_new, order, prefix="a")
