"""
Nonlinear deformation of surface to target space

The purpose of the following script is to deform a list or surface meshs to a target space using
a coordinate mapping.

created by Daniel Haenelt
Date created: 07-02-2019
Last modified: 18-02-2019
"""
import os
from lib.mapping.deform_surface import deform_surface

# input files
input_surf = ["/home/daniel/Schreibtisch/data/sampling_data/layer/lh.layer0",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/lh.layer1",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/lh.layer2",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/lh.layer3",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/lh.layer4",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/lh.layer5",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/lh.layer6",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/lh.layer7",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/lh.layer8",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/lh.layer9",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/rh.layer0",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/rh.layer1",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/rh.layer2",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/rh.layer3",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/rh.layer4",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/rh.layer5",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/rh.layer6",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/rh.layer7",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/rh.layer8",
              "/home/daniel/Schreibtisch/data/sampling_data/layer/rh.layer9",
              ]

input_orig = "/home/daniel/Schreibtisch/data/sampling_data/orig.mgz"
input_deform = "/home/daniel/Schreibtisch/data/sampling_data/ge_epi2/epi2ana.nii.gz"
input_target = "/home/daniel/Schreibtisch/data/sampling_data/meanudata.nii"
path_output = "/home/daniel/Schreibtisch/surf_def"

""" do not edit below """

for i in range(len(input_surf)):
    
    # get hemisphere from input surface file name
    hemi = os.path.splitext(os.path.basename(input_surf[i]))[0]
    
    # deform surface
    deform_surface(input_surf[i], input_orig, input_deform, input_target, hemi, path_output)