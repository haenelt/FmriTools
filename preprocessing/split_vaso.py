"""
Split VASO and BOLD

This scripts splits the vaso and bold time series into separate files. Additionally, volumes at the
beginning (non steady-state volumes) and at the end can be discarded. Note that number of volumes is
the individual volumes, i.e. not the mutual volume TR consisting of both vaso and bold volume.

created by Daniel Haenelt
Date created: 02-05-2018             
Last modified: 02-05-2019  
"""
import os
import numpy as np
import nibabel as nb

# input data
img_input = ["/nobackup/actinium2/haenelt/VasoTest/flicker/Run_1/data.nii",
             "/nobackup/actinium2/haenelt/VasoTest/flicker/Run_2/data.nii",
             "/nobackup/actinium2/haenelt/VasoTest/flicker/Run_3/data.nii",
             "/nobackup/actinium2/haenelt/VasoTest/flicker/Run_4/data.nii",
             "/nobackup/actinium2/haenelt/VasoTest/flicker/Run_5/data.nii",
             "/nobackup/actinium2/haenelt/VasoTest/flicker/Run_6/data.nii",
             ]

# paramteers
start_vol = 10
end_vol = 3

""" do not edit below """

for i in range(len(img_input)):
    
    # load data
    data = nb.load(img_input[i])
    data_array = data.get_fdata()
    
    # discard volumes
    data_array = data_array[:,:,:,start_vol:-end_vol]
    
    # split into even and odd runs
    t_even = np.arange(0,np.shape(data_array)[3],2)
    t_odd = np.arange(1,np.shape(data_array)[3],2)
    
    vaso_array = data_array[:,:,:,t_even]
    bold_array = data_array[:,:,:,t_odd]
    
    # new array length
    data.header["dim"][4] = np.shape(vaso_array)[3]
    
    output = nb.Nifti1Image(vaso_array,data.affine,data.header)
    nb.save(output,os.path.join(os.path.dirname(img_input[i]),"vaso_basis.nii"))
    
    output = nb.Nifti1Image(bold_array,data.affine,data.header)
    nb.save(output,os.path.join(os.path.dirname(img_input[i]),"bold_basis.nii"))