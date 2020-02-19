"""
Bold correction of VASO data

This scripts corrects a vaso time series for bold contamination. First, both time series are 
upsampled to a common time grid. BOLD correction is performed by dividing both time series. In the 
end, unrealistic vaso values are removed.

Before running the script, login to queen via ssh and set the afni environment by calling AFNI in 
the terminal.

created by Daniel Haenelt
Date created: 19-02-2020
Last modified: 19-02-2020
"""
import os
import numpy as np
import nibabel as nb
from lib.io.get_filename import get_filename
from lib.utils.regrid_time_series import regrid_time_series

# input data
img_vaso = ["/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_1/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_2/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_3/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_4/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_5/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_6/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_7/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_8/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_9/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_10/uvaso.nii",
            ]

img_bold = ["/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_1/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_2/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_3/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_4/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_5/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_6/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_7/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_8/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_9/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_10/ubold.nii",
            ]

# parameters
TR = 5 # effective TR of bold+vaso
vaso_shift = 2.8 # start of vaso block (asymmetric TR)
vaso_threshold = 6

""" do not edit below """

for i in range(len(img_vaso)):
    
    # get filenames
    path_bold, name_bold, ext_bold = get_filename(img_bold[i])
    path_vaso, name_vaso, ext_vaso = get_filename(img_vaso[i])
    
    # upsample time series
    regrid_time_series(img_bold[i], path_bold[i], TR, TR/2, t_start=0)    
    regrid_time_series(img_vaso[i], path_vaso[i], TR, TR/2, t_start=vaso_shift)    

    # new filenames
    file_bold = os.path.join(path_bold,name_bold+"_upsampled"+ext_bold)
    file_vaso = os.path.join(path_vaso,name_vaso+"_upsampled"+ext_vaso)
    file_vaso_corrected = os.path.join(path_vaso,name_vaso+"_upsampled_corrected"+ext_vaso)
    
    # load bold data
    bold = nb.load(file_bold)
    bold_array = bold.get_fdata()
    
    # load vaso data
    vaso = nb.load(file_vaso)
    vaso_array = vaso.get_fdata()
    
    # overwrite vaso files at the beginning
    n = 0
    while True:
        if np.sum(vaso_array[:,:,:,n]) == 0:
            n += 1
        else:
            break
    
    for i in range(n):
        vaso_array[:,:,:,i] = vaso_array[:,:,:,n]

    # bold correction
    vaso_array = np.divide(vaso_array, bold_array)

    # clean vaso data that are unrealistic
    vaso_array[vaso_array < 0] = 0
    vaso_array[vaso_array >= vaso_threshold] = vaso_threshold

    # write output
    output = nb.Nifti1Image(vaso_array, vaso.affine, vaso.header)
    nb.save(output, file_vaso_corrected)

    # change TR in header
    os.system("3drefit " + \
              "-TR " + str(TR/2) + " " + \
              file_bold)
    
    os.system("3drefit " + \
              "-TR " + str(TR/2) + " " + \
              file_vaso)

    os.system("3drefit " + \
              "-TR " + str(TR/2) + " " + \
              file_vaso_corrected)