# -*- coding: utf-8 -*-

# external inputs
import numpy as np
import nibabel as nb
from nibabel.freesurfer.io import read_geometry
from nibabel.freesurfer.io import write_geometry
from gbb.normal.get_normal import get_normal


"""
Relief plot on surface mesh

In the following script, a contrast is plotted as a relief on a corresponding
surface mesh. The contrast data is filtered by a sigmoid function.

created by Daniel Haenelt
Date created: 06-03.2020         
Last modified: 12-10-2020  
"""

file_surf = "/home/daniel/source/BlenderCBS-master/Haenelt/data/lh.refined_enhanced_inflated"
file_data = "/home/daniel/source/BlenderCBS-master/Haenelt/data/lh.spmT_left_right_GE_EPI2_upsampled_avg_layer4_16.mgh"
file_out = "/home/daniel/source/BlenderCBS-master/Haenelt/data/lh.refined_enhanced_inflated_relief"
scale_factor = 2

# do not edit below

def sigmoid(x):
  return 1 / (1 + np.exp(-x))

# load geometry
vtx, fac = read_geometry(file_surf)

# load contrast
data = np.zeros_like(vtx)
data[:,0] = nb.load(file_data).get_fdata()[:,0,0]
data[:,1] = nb.load(file_data).get_fdata()[:,0,0]
data[:,2] = nb.load(file_data).get_fdata()[:,0,0]

# get normals
normal = get_normal(vtx, fac)

# make relief
vtx += scale_factor * sigmoid(data) * normal

# write output
write_geometry(file_out, vtx, fac)
