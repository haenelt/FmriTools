# -*- coding: utf-8 -*-
"""
Temporal correction

This scripts performs slice timing correction and temporal regridding on a nifti
time series. The script needs an installation of afni.

"""

# local inputs
from fmri_tools.preprocessing.slice_timing_correction import slice_timing_correction

# input
file_in = [
    "/data/pt_01880/Experiment1_ODC/p4/retinotopy/pol_anticlock/data.nii",
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

for i in range(len(file_in)):
    slice_timing_correction(file_in[i], TR_old, TR_new, order, prefix="a")
