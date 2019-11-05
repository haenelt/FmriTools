
"""
Split VASO and BOLD

This scripts splits the vaso and bold time series into separate files. The time series always starts
with the non-nulled BOLD volume. So, all even and odd time points are denoted as bold and vaso,
respectively. Additionally, volumes at the beginning (non steady-state volumes) are overwritten by 
following steady-state volumes. Furthermore, volumes at the end of the scan can be discarded. Note 
that the number of volumes refers here to the individual volumes in the bold+vaso time series.

created by Daniel Haenelt
Date created: 02-05-2018             
Last modified: 20-10-2019  
"""
import os
import numpy as np
import nibabel as nb

# input data
img_input = ["/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_1/data.nii",
             "/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_2/data.nii",
             "/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_3/data.nii",
             "/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_4/data.nii",
             "/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_5/data.nii",
             "/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_6/data.nii",
             "/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_7/data.nii",
             "/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_8/data.nii",
             "/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_9/data.nii",
             "/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_10/data.nii",
             ]

# paramteers
start_vol = 2
end_vol = 5

""" do not edit below """

for i in range(len(img_input)):
    
    # load data
    data = nb.load(img_input[i])
    data_array = data.get_fdata()
    
    # overwrite non steady-state volumes
    data_array[:,:,:,0:start_vol] = data_array[:,:,:,start_vol:2*start_vol]
    
    # discard volumes at the end
    if end_vol != 0:
        data_array = data_array[:,:,:,:-end_vol]
    
    # split into even and odd runs
    t_even = np.arange(0,np.shape(data_array)[3],2)
    t_odd = np.arange(1,np.shape(data_array)[3],2)
    
    bold_array = data_array[:,:,:,t_even]
    vaso_array = data_array[:,:,:,t_odd]
    
    # new array length
    data.header["dim"][4] = np.shape(vaso_array)[3]
    
    output = nb.Nifti1Image(vaso_array,data.affine,data.header)
    nb.save(output,os.path.join(os.path.dirname(img_input[i]),"vaso.nii"))
    
    output = nb.Nifti1Image(bold_array,data.affine,data.header)
    nb.save(output,os.path.join(os.path.dirname(img_input[i]),"bold.nii"))