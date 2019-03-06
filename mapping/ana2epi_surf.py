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
input_surf = ["/data/pt_01880/V2STRIPES/p7/anatomy/layer/lh.layer0",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/lh.layer1",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/lh.layer2",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/lh.layer3",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/lh.layer4",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/lh.layer5",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/lh.layer6",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/lh.layer7",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/lh.layer8",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/lh.layer9",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/rh.layer0",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/rh.layer1",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/rh.layer2",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/rh.layer3",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/rh.layer4",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/rh.layer5",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/rh.layer6",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/rh.layer7",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/rh.layer8",
              "/data/pt_01880/V2STRIPES/p7/anatomy/layer/rh.layer9",
              ]

input_orig = "/data/pt_01880/V2STRIPES/p7/anatomy/freesurfer/mri/orig.mgz"
input_deform = "/data/pt_01880/V2STRIPES/p7/deformation/colour/ge_epi1/epi2orig.nii.gz"
input_target = "/data/pt_01880/V2STRIPES/p7/colour/GE_EPI1/Run_1/meanudata.nii"
path_output = "/data/pt_01880/V2STRIPES/p7/anatomy/layer_def/colour/ge_epi1"

""" do not edit below """

for i in range(len(input_surf)):
    
    # get hemisphere from input surface file name
    hemi = os.path.splitext(os.path.basename(input_surf[i]))[0]
    
    # deform surface
    deform_surface(input_surf[i], input_orig, input_deform, input_target, hemi, path_output)
