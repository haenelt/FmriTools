# -*- coding: utf-8 -*-
"""
Apply nonlinear deformation to surface mesh

The purpose of the following script is to apply two deformations in succession
(orig -> ana -> epi) to a surface mesh using generated coordinate mappings. The
script needs an installation of freesurfer.

"""

import os
from os.path import basename, join

from ..registration.surf import deform_surface

# input files
input_surf = ["/data/pt_01880/temp_reg/data/lh.layer10"]
input_orig = "/data/pt_01880/temp_reg/data/orig.mgz"
input_ana = "/data/pt_01880/temp_reg/data/S22_MP2RAGE_0p7_T1_Images_2.45.nii"
input_epi = "/data/pt_01880/temp_reg/data/mean_test2.nii"
input_mask = ""  # epi2ana mask
input_deform1 = "/data/pt_01880/temp_reg/deformation/header/S22_MP2RAGE_0p7_T1_Images_2.45_2_orig_scanner.nii"  # ana2orig
input_deform2 = "/data/pt_01880/temp_reg/deformation/enhanced/epi2ana.nii.gz"  # epi2ana
path_output = "/data/pt_01880/temp_reg/surf/enhanced"

# do not edit below

for i in range(len(input_surf)):
    # orig -> ana
    deform_surface(
        input_surf[i],
        input_orig,
        input_deform1,
        input_ana,
        path_output,
        input_mask=None,
        interp_method="trilinear",
        smooth_iter=0,
        flip_faces=True,
        cleanup=False,
    )

    # rename output
    os.rename(
        join(path_output, basename(input_surf[i]) + "_def"),
        join(path_output, basename(input_surf[i]) + "_def1"),
    )

    # ana -> epi
    deform_surface(
        join(path_output, basename(input_surf[i]) + "_def1"),
        input_ana,
        input_deform2,
        input_epi,
        path_output,
        input_mask,
        interp_method="trilinear",
        smooth_iter=0,
        flip_faces=False,
        cleanup=False,
    )

    # rename output
    os.rename(
        join(path_output, basename(input_surf[i]) + "_def1_def"),
        join(path_output, basename(input_surf[i]) + "_def2"),
    )
    os.rename(
        join(path_output, basename(input_surf[i]) + "_def1_ind.txt"),
        join(path_output, basename(input_surf[i]) + "_def2_ind"),
    )
