# -*- coding: utf-8 -*-

# python standard library inputs
# external inputs
# local inputs


"""
bla...

# match vertex numbers
# mris_extract main_component
# inflated
# meshlines
# layers

created by Daniel Haenelt
Date created: 24-10-2020
Last modified: 24-10-2020
"""
import os
import sys
import subprocess
import numpy as np
from nibabel.freesurfer.io import read_geometry
from fmri_tools.surface import match_vertex_number
from fmri_tools.surface import inflate_surf_mesh
from fmri_tools.io import get_filename
from nibabel.freesurfer.io import write_geometry

input_pial_surf = "/data/pt_01880/Experiment2_Rivalry/p1/anatomy/gbb/lh.pial_def2_refined"
input_white_surf = "/data/pt_01880/Experiment2_Rivalry/p1/anatomy/gbb/lh.white_def2_refined"
input_pial_ind = "/data/pt_01880/Experiment2_Rivalry/p1/anatomy/dense_epi/lh.pial_def2_ind"
input_white_ind = "/data/pt_01880/Experiment2_Rivalry/p1/anatomy/dense_epi/lh.white_def2_ind"
input_pial_ind2 = "/data/pt_01880/Experiment2_Rivalry/p1/anatomy/gbb/lh.pial_def2_refined_ind.txt"
input_white_ind2 = "/data/pt_01880/Experiment2_Rivalry/p1/anatomy/gbb/lh.white_def2_refined_ind.txt"

path_output = ""
niter_inflate = 30

# do not edit below

hemi = ["lh","rh"]
i = 0

# load surfaces
vtx_white, fac_white = read_geometry(input_white_surf)
vtx_pial, _ = read_geometry(input_pial_surf)

# load index files
ind_white = np.loadtxt(input_white_ind).astype(int)
ind_pial = np.loadtxt(input_pial_ind).astype(int)
ind_white2 = np.loadtxt(input_white_ind2).astype(int)
ind_pial2 =np.loadtxt(input_pial_ind2).astype(int)

ind_white = ind_white[ind_white2]
ind_pial = ind_pial[ind_pial2]

# get filenames
_, _, name_white = get_filename(input_white_surf)
_, _, name_pial = get_filename(input_pial_surf)

# make folders
path_dense = os.path.join(path_output, "dense_refined")
path_layer = os.path.join(path_output, "layer")

if not os.path.exists(path_output):
    os.makedirs(path_output)

if not os.path.exists(path_dense):
    os.mkdir(path_dense)
else:
    raise FileExistsError("Dense folder already exists!")

if not os.path.exists(path_layer):
    os.makedirs(path_layer)
else:
    raise FileExistsError("Layer folder already exists!")

#%% match vertex numbers
vtx_white, vtx_pial, fac_white, ind = match_vertex_number(vtx_white, 
                                                          vtx_pial, 
                                                          fac_white, 
                                                          ind_white, 
                                                          ind_pial)

file_white = os.path.join(path_dense, hemi[i]+"."+name_white)
write_geometry(file_white, vtx_white, fac_white)

file_pial = os.path.join(path_dense, hemi[i]+"."+name_pial)
write_geometry(file_pial, vtx_pial, fac_white)

#%% extract main component
try:
    subprocess.run(['mris_extract_main_component', 
                    file_white, 
                    file_white+"_main"], check = True)
except subprocess.CalledProcessError:
    sys.exit("Main component extraction failed!")

try:
    subprocess.run(['mris_extract_main_component', 
                    file_pial, 
                    file_pial+"_main"], check = True)
except subprocess.CalledProcessError:
    sys.exit("Main component extraction failed!")

#%% inflation
inflate_surf_mesh(file_white+"_main", 
                  file_white+"_inflate"+str(niter_inflate), 
                  niter_inflate)

