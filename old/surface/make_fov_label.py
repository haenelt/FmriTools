# -*- coding: utf-8 -*-
"""
Create FOV ROI

This script computes a label file to include only vertices (in further analysis) which
are within the intersection of multiple volumes.

"""

import nibabel as nb
import numpy as np

from ..io.surf import write_label
from ..surface.label import roi_fov
from ..surface.mesh import Mesh

# input
file_surf = [
    "/data/pt_01880/Experiment1_ODC/p1/anatomy/layer/lh.layer_0",
    "/data/pt_01880/Experiment1_ODC/p1/anatomy/layer/lh.layer_1",
    "/data/pt_01880/Experiment1_ODC/p1/anatomy/layer/lh.layer_2",
    "/data/pt_01880/Experiment1_ODC/p1/anatomy/layer/lh.layer_3",
    "/data/pt_01880/Experiment1_ODC/p1/anatomy/layer/lh.layer_4",
    "/data/pt_01880/Experiment1_ODC/p1/anatomy/layer/lh.layer_5",
    "/data/pt_01880/Experiment1_ODC/p1/anatomy/layer/lh.layer_6",
    "/data/pt_01880/Experiment1_ODC/p1/anatomy/layer/lh.layer_7",
    "/data/pt_01880/Experiment1_ODC/p1/anatomy/layer/lh.layer_8",
    "/data/pt_01880/Experiment1_ODC/p1/anatomy/layer/lh.layer_9",
    "/data/pt_01880/Experiment1_ODC/p1/anatomy/layer/lh.layer_10",
]
file_cmap = [
    "/data/pt_01880/Experiment1_ODC/p1/deformation/odc/GE_EPI3/source2target.nii.gz",
    "/data/pt_01880/Experiment1_ODC/p1/deformation/odc/GE_EPI4/source2target.nii.gz",
    "/data/pt_01880/Experiment1_ODC/p1/deformation/odc/SE_EPI1/source2target.nii.gz",
    "/data/pt_01880/Experiment1_ODC/p1/deformation/odc/SE_EPI2/source2target.nii.gz",
    "/data/pt_01880/Experiment1_ODC/p1/deformation/odc/VASO1/source2target.nii.gz",
    "/data/pt_01880/Experiment1_ODC/p1/deformation/odc/VASO2/source2target.nii.gz",
]
file_target = [
    "/data/pt_01880/Experiment1_ODC/p1/odc/GE_EPI3/diagnosis/mean_udata.nii",
    "/data/pt_01880/Experiment1_ODC/p1/odc/GE_EPI4/diagnosis/mean_udata.nii",
    "/data/pt_01880/Experiment1_ODC/p1/odc/SE_EPI1/diagnosis/mean_udata.nii",
    "/data/pt_01880/Experiment1_ODC/p1/odc/SE_EPI2/diagnosis/mean_udata.nii",
    "/data/pt_01880/Experiment1_ODC/p1/odc/VASO1/diagnosis/mean_ubold.nii",
    "/data/pt_01880/Experiment1_ODC/p1/odc/VASO2/diagnosis/mean_ubold.nii",
]
file_out = "/data/pt_01880/zzz.label"

# do not edit below

# initialize label array
arr_label = np.arange(len(Mesh.from_file(file_surf[0]).verts))

# loop through target files and layers
for ft, fc in zip(file_target, file_cmap):
    hdr = nb.load(ft).header
    vol_dims = hdr["dim"][1:4]
    vol_ds = hdr["pixdim"][1:4]
    for fs in file_surf:
        mesh = Mesh.from_file(fs)
        # transform vertices to voxel target space
        mesh.transform_coords(fc, ft)
        # new label from intersection with old label
        arr_label = np.intersect1d(arr_label, roi_fov(mesh.verts, vol_dims, vol_ds))

# write output
write_label(file_out, arr_label)
