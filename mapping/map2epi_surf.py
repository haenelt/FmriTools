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
input_surf = "/home/daniel/Schreibtisch/test/lh.layer0_def"
input_vol = "/home/daniel/Schreibtisch/surface_voxel.nii"
path_output = "/home/daniel/Schreibtisch"
input_white = "/home/daniel/Schreibtisch/data/sampling_data/lh.white"
input_ind = "/home/daniel/Schreibtisch/test/lh.layer0_ind.txt"
cleanup=False

""" do not edit below """

# get hemisphere from input surface file name
hemi = os.path.splitext(os.path.basename(input_surf))[0]
    
# deform surface
map2surface(input_surf, input_vol, hemi, path_output, input_white, input_ind, cleanup)