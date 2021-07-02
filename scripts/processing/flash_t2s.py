# -*- coding: utf-8 -*-
"""
Flash T2s fitting

The script applies the nighres t2s fitting module to a data set.

"""

# python standard library inputs
import os

# external inputs
import nibabel as nb
from nighres.intensity import flash_t2s_fitting

# local inputs
from fmri_tools.io.get_filename import get_filename

# input
file_list = [
    "/data/pt_01880/Experiment1_ODC/p4/anatomy/flash3/S13_3D_GRE_3ech_iso0p5_slab_8.42.nii",
    "/data/pt_01880/Experiment1_ODC/p4/anatomy/flash3/S13_3D_GRE_3ech_iso0p5_slab_16.03.nii",
    "/data/pt_01880/Experiment1_ODC/p4/anatomy/flash3/S13_3D_GRE_3ech_iso0p5_slab_25.nii"]
te_list = [8.42, 16.03, 25.00]  # in ms
name_output = "3D_GRE_3ech_iso0p5_slab"

# do not edit below

# get output path from first input entry
path_output, _, _ = get_filename(file_list[0])

# t2s fitting
res = flash_t2s_fitting(file_list, te_list)

# write output
nb.save(res["t2s"], os.path.join(path_output, name_output + "_t2s.nii"))
nb.save(res["r2s"], os.path.join(path_output, name_output + "_r2s.nii"))
nb.save(res["s0"], os.path.join(path_output, name_output + "_s0.nii"))
nb.save(res["residuals"], os.path.join(path_output,
                                       name_output + "_residuals.nii"))
