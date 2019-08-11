"""
Split VASO and BOLD

This scripts splits the vaso and bold time series into separate files. Additionally, volumes at the
beginning (non steady-state volumes) are overwritten by following steady-state volumes. Furthermore,
volumes at the end of the scan can be discarded. Note that the number of volumes refers here to the
individual volumes in the bold+vaso time series.

created by Daniel Haenelt
Date created: 02-05-2018             
Last modified: 02-05-2019  
"""
import os
import numpy as np
import nibabel as nb

# input data
img_input = ["/nobackup/actinium2/haenelt/phantom_20190809/Run_2/data.nii",
             ]

# paramteers
start_vol = 10
end_vol = 3

""" do not edit below """

for i in range(len(img_input)):
    
    # load data
    data = nb.load(img_input[i])
    data_array = data.get_fdata()
    
    # overwrite non steady-state volumes
    data_array[:,:,:,0:start_vol] = data_array[:,:,:,start_vol:2*start_vol]
    
    # discard volumes at the end
    data_array = data_array[:,:,:,:-end_vol]
    
    # split into even and odd runs
    t_even = np.arange(0,np.shape(data_array)[3],2)
    t_odd = np.arange(1,np.shape(data_array)[3],2)
    
    vaso_array = data_array[:,:,:,t_even]
    bold_array = data_array[:,:,:,t_odd]
    
    # new array length
    data.header["dim"][4] = np.shape(vaso_array)[3]
    
    output = nb.Nifti1Image(vaso_array,data.affine,data.header)
    nb.save(output,os.path.join(os.path.dirname(img_input[i]),"vaso.nii"))
    
    output = nb.Nifti1Image(bold_array,data.affine,data.header)
    nb.save(output,os.path.join(os.path.dirname(img_input[i]),"bold.nii"))