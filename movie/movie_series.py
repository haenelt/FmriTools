"""
Make time series movie

The purpose of the following script is to make a movie to illustrate directly the demeaned time 
series on the cortical surface.

created by Daniel Haenelt
Date created: 12-11-2019
Last modified: 12-11-2019
"""
import os
import shutil as sh
import numpy as np
import nibabel as nb
from nighres.registration import apply_coordinate_mappings
from lib.utils.get_mean4d import get_mean4d
from lib.processing.demean_time_series import demean_time_series
from lib.mapping import map2surface

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

# path to SPM12 folder
pathSPM = "/data/pt_01880/source/spm12"
pathLIB = "/data/hu_haenelt/projects/scripts/lib/preprocessing"

""" do not edit below """

# prefix
tmp = np.random.randint(0, 10, 5)
prefix = ''.join(str(i) for i in tmp)+'_'

# change to lib folder
os.chdir(pathLIB)

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
    os.system("matlab" + \
              " -nodisplay -nodesktop -r " + \
              "\"baseline_correction(\'{0}\', {1}, {2}, \'{3}\', \'{4}\'); exit;\"". \
              format(input_series[i], TR, cutoff_highpass, pathSPM, prefix))
    
    # move baseline corrected time series to output folder
    sh.move(os.path.join(os.path.dirname(input_series[i]),prefix+os.path.basename(input_series[i])),
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
        
        # hemisphere
        hemi = os.path.splitext(os.path.basename(input_surf[j]))[0]
    
        # sample on surface
        map2surface(input_surf[j], 
                    os.path.join(path_def,str(i+1)+"_def-img.nii.gz"), 
                    hemi, 
                    path_surf, 
                    input_white=None, 
                    input_ind=None, 
                    cleanup=True)
