# -*- coding: utf-8 -*-
"""
Percent signal change conversion

This scripts converts timeseries data into percent signal change.

"""

# local inputs
from fmri_tools.preprocessing.scale_timeseries import scale_timeseries

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
cutoff_psc = 50

# do not edit below

for i in range(len(file_in)):
    scale_timeseries(file_in[i], cutoff=cutoff_psc, prefix="p")
