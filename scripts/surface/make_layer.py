# -*- coding: utf-8 -*-

# python standard library inputs
# external inputs
# local inputs


"""
bla...

created by Daniel Haenelt
Date created: xx-xx-xxxx
Last modified: xx-xx-xxxx
"""
import os
import numpy as np
from nibabel.freesurfer.io import read_geometry
from fmri_tools.surface import match_vertex_number
from fmri_tools.surface import inflate_surf_mesh
from nibabel.freesurfer.io import write_geometry

# inflate
# match vertex numbers
# mris_extract_main_component
# meshlines
# layers

input_pial_surf = "/data/pt_01880/Experiment2_Rivalry/p1/anatomy/gbb/lh.pial_def2_refined"
input_white_surf = "/data/pt_01880/Experiment2_Rivalry/p1/anatomy/gbb/lh.white_def2_refined"
input_pial_ind = "/data/pt_01880/Experiment2_Rivalry/p1/anatomy/dense_epi/lh.pial_def2_ind"
input_white_ind = "/data/pt_01880/Experiment2_Rivalry/p1/anatomy/dense_epi/lh.white_def2_ind"
input_pial_ind2 = "/data/pt_01880/Experiment2_Rivalry/p1/anatomy/gbb/lh.pial_def2_refined_ind.txt"
input_white_ind2 = "/data/pt_01880/Experiment2_Rivalry/p1/anatomy/gbb/lh.white_def2_refined_ind.txt"

vtx_white, fac = read_geometry(input_white_surf)
vtx_pial, fac2 = read_geometry(input_pial_surf)

ind_white = np.loadtxt(input_white_ind).astype(int)
ind_pial = np.loadtxt(input_pial_ind).astype(int)
ind_white2 = np.loadtxt(input_white_ind2).astype(int)
ind_pial2 =np.loadtxt(input_pial_ind2).astype(int)

ind_white = ind_white[ind_white2]
ind_pial = ind_pial[ind_pial2]

#%%

v_white, v_pial, f, ind = match_vertex_number(vtx_white, 
                                              vtx_pial, 
                                              fac, 
                                              ind_white, 
                                              ind_pial)


#%%

write_geometry("/data/pt_01880/test", v_white, f)
write_geometry("/data/pt_01880/test2", v_pial, f)

#%%

file_in = "/data/pt_01880/test_main"


inflate_surf_mesh(file_in, file_in+"_ind", 30)

