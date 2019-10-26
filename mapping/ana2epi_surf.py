"""
Nonlinear deformation of surface to target space

The purpose of the following script is to deform a list or surface meshs to a target space using
a coordinate mapping.

Before running the script, login to queen via ssh and set the freesurfer environment by calling 
FREESURFER in the terminal.

created by Daniel Haenelt
Date created: 07-02-2019
Last modified: 18-02-2019
"""
import os
from lib.mapping.deform_surface import deform_surface

# input files
input_surf = ["/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/lh.layer5",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/rh.layer5",
              ]

input_orig = "/data/pt_01880/Experiment1_ODC/p3/anatomy/freesurfer/mri/orig.mgz"
input_deform = "/data/pt_01880/Experiment1_ODC/p3/deformation/odc/GE_EPI1/epi2orig.nii.gz"
input_target = "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI1/diagnosis/mean_data.nii"
path_output = "/data/pt_01880/test/non_rigid"

""" do not edit below """

for i in range(len(input_surf)):
    
    # get hemisphere from input surface file name
    hemi = os.path.splitext(os.path.basename(input_surf[i]))[0]
    
    # deform surface
    deform_surface(input_surf[i], input_orig, input_deform, input_target, hemi, path_output)
