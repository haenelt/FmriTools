# -*- coding: utf-8 -*-

# python standard library inputs
import os
import shutil as sh
from os.path import join, exists

# external inputs
from scipy.io import loadmat, savemat
from nibabel.affines import apply_affine
from nibabel.freesurfer.io import write_geometry
from gbb.utils.vox2ras import vox2ras


"""
Recursive Boundary-Based Registration

This script performs recurive BBR to register a functional data set with 
surfaces generated from a separate whole-brain anatomy.

parameter mode (explanation from one of Tim vanMourik's scripts):
`r' for rotation around all axes, `rx', `ry', or `rz' for a single rotation 
around the respective axis and for example `rxrz' for a rotation around the 
x-axis and the z-axis. In the same way the scaling and translation can be 
modified. For example the combination `rystytz' would optimise for a rotation 
around the y-axis, a scaling in all direction and a translation along the y-axis 
and along the z-axis. The default is 'rst' for all transformations.

The script needs an installation of freesurfer.

created by Daniel Haenelt
Date created: 12-12-2019             
Last modified: 13-10-2020  
"""

# input surface
lh_white = "/data/pt_01880/odc_temp/surface/match/lh.layer10_def2_smooth_match"
rh_white = "/data/pt_01880/odc_temp/surface/match/rh.layer10_def2_smooth_match"
lh_pial = "/data/pt_01880/odc_temp/surface/match/lh.layer0_def2_smooth_match"
rh_pial = "/data/pt_01880/odc_temp/surface/match/rh.layer0_def2_smooth_match"

# parameters
subjects_dir = "/data/pt_01880/odc_temp/surface/rBBR" # absolute path to subject directory
ref_vol = "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/diagnosis/mean_data.nii" # absolute path to reference volume
mask = "" # absolute path to mask volume
min_vox = 10 
min_vtx = 500
cuboid = True
tetrahedra = False
nn_smooth = 0.5
registration_mode = "syty"
reverse_contrast = True # if inner surface is darker than outer surface
cost_method = "GreveFischl" # GreveFischl, sum, sumOfSquares
contrast_method = "gradient" # gradient, average, fixedDistance, relativeDistance, extrema

# add spm, fmri_tools and OpenFmriAnalysis to path
pathSPM = "/data/pt_01880/source/spm12"
pathMOURIK = "/data/hu_haenelt/projects/OpenFmriAnalysis"
pathFMRITOOLS = "/data/hu_haenelt/projects/FmriTools/fmri_tools"

# do not edit below

# make output folder
if not exists(subjects_dir):
    os.makedirs(subjects_dir)

path_input = join(subjects_dir, "input")
if not exists(path_input):
    os.makedirs(path_input)

# copy input files to input folder
sh.copyfile(lh_white, join(path_input,"lh.white"))
sh.copyfile(rh_white, join(path_input,"rh.white"))
sh.copyfile(lh_pial, join(path_input,"lh.pial"))
sh.copyfile(rh_pial, join(path_input,"rh.pial"))
sh.copyfile(ref_vol, join(path_input,"reference_volume.nii"))
sh.copyfile(mask, join(path_input,"mask.nii")) if len(mask) else None

# filenames
input_surf_ras = "boundaries_in_ras.mat"
input_surf_vox = "boundaries_in_vox.mat"
output_surf_ras = "boundaries_out_ras.mat"
output_surf_vox = "boundaries_out_vox.mat"
output_cmap = "cmap_rBBR.nii"

lh_white = join(path_input,"lh.white")
rh_white = join(path_input,"rh.white")
lh_pial = join(path_input,"lh.pial")
rh_pial = join(path_input,"rh.pial")
ref_vol = join(path_input,"reference_volume.nii")
in_surf_mat_ras = join(path_input, input_surf_ras)
in_surf_mat_vox = join(path_input, input_surf_vox)
out_surf_mat_ras = join(subjects_dir, output_surf_ras)
out_surf_mat_vox = join(subjects_dir, output_surf_vox)
input_cfg = join(path_input,"cfg.mat")

# get affine vox2ras-tkr and ras2vox-tkr transformation to reference volume
vox2ras_tkr, ras2vox_tkr = vox2ras(ref_vol)

# surf2mat
cwd = os.getcwd()
os.chdir(join(pathFMRITOOLS,"io"))
os.system("matlab" + \
          " -nodisplay -nodesktop -r " + \
          "\"surf2mat(\'{0}\', \'{1}\', \'{2}\', \'{3}\', \'{4}\'); exit;\"". \
          format(lh_white, rh_white, lh_pial, rh_pial, in_surf_mat_ras))
os.chdir(cwd)
    
# ras2vox
data = loadmat(in_surf_mat_ras)

# apply ras2vox transformation to vertices
for i in range(2):
    data["wSurface"][0][i] = apply_affine(ras2vox_tkr, data["wSurface"][0][i])
    data["pSurface"][0][i] = apply_affine(ras2vox_tkr, data["pSurface"][0][i])
    
# save surfaces in voxel space
savemat(in_surf_mat_vox, data)

# get configuration mat-file for rBBR
cfg = {}
cfg["subjects_dir"] = subjects_dir
cfg["reference_volume"] = join("input","reference_volume.nii")
cfg["input_surf"] = join("input",input_surf_vox)
cfg["mask"] = join("input","mask.nii") if len(mask) else ""
cfg["min_vox"] = float(min_vox)
cfg["min_vtx"] = float(min_vtx) 
cfg["cuboid"] = cuboid
cfg["tetrahedra"] = tetrahedra 
cfg["nn_smooth"] = float(nn_smooth)
cfg["registration_mode"] = registration_mode
cfg["reverse_contrast"] = reverse_contrast
cfg["cost_method"] = cost_method
cfg["contrast_method"] = contrast_method
cfg["output_surf"] = output_surf_vox
cfg["output_cmap"] = output_cmap

savemat(input_cfg, cfg)

# run rBBR
cwd = os.getcwd()
os.chdir(join(pathFMRITOOLS,"registration"))
os.system("matlab" + \
          " -nodisplay -nodesktop -r " + \
          "\"run_rBBR(\'{0}\', \'{1}\', \'{2}\'); exit;\"". \
          format(input_cfg, pathSPM, pathMOURIK))
os.chdir(cwd)

# vox2ras
data = loadmat(out_surf_mat_vox)

# apply ras2vox transformation to output
for i in range(2):
    data["wSurface"][0][i] = data["wSurface"][0][i][:,0:3]
    data["pSurface"][0][i] = data["pSurface"][0][i][:,0:3]  
    
    data["wSurface"][0][i] = apply_affine(vox2ras_tkr, data["wSurface"][0][i])
    data["pSurface"][0][i] = apply_affine(vox2ras_tkr, data["pSurface"][0][i])
    
# save matfile
savemat(out_surf_mat_ras, data)

# save surfaces in freesurfer format
write_geometry(join(subjects_dir,"lh.white"), data["wSurface"][0][0], data["faceData"][0][0])
write_geometry(join(subjects_dir,"rh.white"), data["wSurface"][0][1], data["faceData"][0][1])
write_geometry(join(subjects_dir,"lh.pial"), data["pSurface"][0][0], data["faceData"][0][0])
write_geometry(join(subjects_dir,"rh.pial"), data["pSurface"][0][1], data["faceData"][0][1])
