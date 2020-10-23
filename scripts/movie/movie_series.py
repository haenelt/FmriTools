# -*- coding: utf-8 -*-

# python standard library inputs
import os
import shutil as sh

# external inputs
import numpy as np
import nibabel as nb
from nighres.registration import apply_coordinate_mappings

# local inputs
from fmri_tools.utils.get_mean4d import get_mean4d
from fmri_tools.processing.demean_time_series import demean_time_series
from fmri_tools.mapping import map2surface


"""
Make time series movie

The purpose of the following script is to make a movie to illustrate directly 
the demeaned time series on the cortical surface.

created by Daniel Haenelt
Date created: 12-11-2019
Last modified: 23-10-2020
"""

input_series = ["/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_1/udata.nii",
                "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_2/udata.nii",
                "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_3/udata.nii",
                "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_4/udata.nii",
                "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_5/udata.nii",
                "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_6/udata.nii",
                "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_7/udata.nii",
                "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_8/udata.nii",
                "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_9/udata.nii",
                "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_10/udata.nii",
                ]
input_deformation = "/data/pt_01880/Experiment1_ODC/p3/deformation/odc/GE_EPI2_rigid/epi2orig.nii.gz"
input_surf = ["/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/lh.layer5",
              ]
path_output = "/data/pt_01880/odc_movie"

# paramteers
TR = 3
cutoff_highpass = 270
start_vol = 0
end_vol = 0
interpolation = "nearest" # can be linear or nearest

# add spm and fmri_tools to path
pathSPM = "/data/pt_01880/source/spm12"
pathFMRITOOLS = "/data/hu_haenelt/projects/FmriTools/fmri_tools"

# do not edit below

# prefix
tmp = np.random.randint(0, 10, 5)
prefix = ''.join(str(i) for i in tmp)+'_'

# change to preprocessing folder in fmri_tools
os.chdir(os.path.join(pathFMRITOOLS,"preprocessing"))

# make output folders
path_series = os.path.join(path_output,"series")
path_native = os.path.join(path_output,"native")
path_def = os.path.join(path_output,"def")
path_surf = os.path.join(path_output,"surf")

if not os.path.exists(path_output):
    os.mkdir(path_output)

if not os.path.exists(path_series):
    os.mkdir(path_series)

if not os.path.exists(path_native):
    os.mkdir(path_native)
    
if not os.path.exists(path_def):
    os.mkdir(path_def)

if not os.path.exists(path_surf):
    os.mkdir(path_surf)

# look for baseline corrected time series
for i in range(len(input_series)):

    path_tmp = os.path.dirname(input_series[i])
    basename_tmp = prefix+os.path.basename(input_series[i])
    if os.path.exists(os.path.join(path_tmp, basename_tmp)):
        raise FileExistsError("Temporary file already exists!")

    os.system("matlab" + \
              " -nodisplay -nodesktop -r " + \
              "\"baseline_correction(\'{0}\', {1}, {2}, \'{3}\', \'{4}\'); exit;\"". \
              format(input_series[i], TR, cutoff_highpass, pathSPM, prefix))
    
    # move baseline corrected time series to output folder
    sh.move(os.path.join(path_tmp, basename_tmp),
            os.path.join(path_series,str(i+1)+".nii"))
    
    # rename input images
    input_series[i] = os.path.join(path_series,str(i+1)+".nii")

# get demeaned mean time series
data = get_mean4d(input_series, path_output="", name_output="", write_output=False)
data = demean_time_series(data, write_output=False)
data_array = data.get_fdata()
    
# discard volumes at the beginning and at the end
if start_vol != 0:
   data_array = data_array[:,:,:,start_vol:]
   
if end_vol != 0:
    data_array = data_array[:,:,:,:-end_vol]

# write single volumes
data.header["dim"][4] = 1
data.header["dim"][0] = 3
for i in range(np.shape(data_array)[3]):
    output = nb.Nifti1Image(data_array[:,:,:,i], data.affine, data.header)
    nb.save(output, os.path.join(path_native,str(i+1)+".nii")) 

# apply deformation    
for i in range(np.shape(data_array)[3]):
    apply_coordinate_mappings(os.path.join(path_native,str(i+1)+".nii"), 
                              input_deformation, 
                              interpolation=interpolation, 
                              padding='closest', 
                              save_data=True, 
                              overwrite=True, 
                              output_dir=path_def,
                              file_name=None,
                              )

# map to ana
for i in range(np.shape(data_array)[3]):
    for j in range(len(input_surf)):
    
        map2surface(input_surf[j], 
                    os.path.join(path_def,str(i+1)+"_def-img.nii.gz"), 
                    True, 
                    path_surf, 
                    interp_method="nearest",
                    input_surf_target=None, 
                    input_ind=None, 
                    cleanup=True)
