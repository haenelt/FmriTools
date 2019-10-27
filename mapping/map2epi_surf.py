"""
Sample volume data on surface

The purpose of the following script is to sample data to a surface in native epi space and map those
data to the surface in conformed freesurfer space.

created by Daniel Haenelt
Date created: 07-02-2019
Last modified: 18-02-2019
"""
import os
from lib.mapping.map2surface import map2surface

# input files
input_surf = "/data/pt_01880/test/fieldmap_and_rigid/rh.layer5_def"
input_vol = "/data/pt_01880/test/fieldmap_and_rigid/results/spmT/native/spmT_left_right_GE_EPI1_fieldmap.nii"
path_output = "/data/pt_01880/test/fieldmap_and_rigid/res"
input_white = "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/rh.layer10"
input_ind = "/data/pt_01880/test/fieldmap_and_rigid/rh.layer5_ind.txt"
cleanup=False

""" do not edit below """

# get hemisphere from input surface file name
hemi = os.path.splitext(os.path.basename(input_surf))[0]
    
# deform surface
map2surface(input_surf, input_vol, hemi, path_output, input_white, input_ind, cleanup)