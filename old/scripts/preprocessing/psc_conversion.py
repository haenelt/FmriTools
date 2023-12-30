# -*- coding: utf-8 -*-
"""
Percent signal change conversion

This scripts converts timeseries data into percent signal change.

"""

import nibabel as nb

from ..io.filename import get_filename
from ..preprocessing.timeseries import ScaleTimeseries

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

for f in file_in:
    data = nb.load(f)
    scale = ScaleTimeseries(data.get_fdata())
    arr_scaled = scale.psc(cutoff_psc)
    output = nb.Nifti1Image(arr_scaled, data.affine, data.header)
    path, basename, ext = get_filename(f)
    nb.save(output, f"{path}/p{basename}{ext}")
