# -*- coding: utf-8 -*-

# python standard library inputs
import os
import shutil as sh

# external inputs
import numpy as np
import nibabel as nb
from sh import gunzip
from nighres.registration import apply_coordinate_mappings


"""
Transform time series to target space

In the following script, epi time series in native space are transformed to a 
target space using a deformation field. The transformed time series get the 
prefix r. Because of heap size limits, the time series is split and the 
deformation is applied separately to each volume.

created by Daniel Haenelt
Date created: 07-08-2019            
Last modified: 23-10-2020  
"""

# input
input_epi = [
        "/data/pt_01880/preprocessing_test/method2/Run_2/udata.nii",
        ]

input_reg = [
        "/data/pt_01880/preprocessing_test/deformation/method2/Run_2/source2target.nii.gz",
        ]

# parameters
interpolation = "linear"
padding = "closest"

# do not edit below

# apply deformation    
if len(input_epi) == len(input_reg):
    for i in range(len(input_epi)):
        
        # make temporary output folder
        tmp = np.random.randint(0, 10, 5)
        tmp_string = ''.join(str(i) for i in tmp)
        path_tmp = os.path.join(os.path.dirname(input_epi[i]),"tmp_"+tmp_string)
        
        if not os.path.exists(path_tmp):
            os.mkdir(path_tmp)
        else:
            raise FileExistsError("Temporary folder already exists!")
        
        # get length of time series
        data = nb.load(input_epi[i])
        data_array = data.get_fdata()
        nt = nb.load(input_epi[i]).header["dim"][4]
        data.header["dim"][4] = 1 # change header for single 3d volumes
        
        for j in range(nt):
            
            # save single time frame
            output = nb.Nifti1Image(data_array[:,:,:,j], data.affine, data.header)
            nb.save(output, os.path.join(path_tmp,str(j)+".nii"))
            
            apply_coordinate_mappings(os.path.join(path_tmp,str(j)+".nii"), 
                                      input_reg[i], 
                                      interpolation=interpolation, 
                                      padding=padding,
                                      save_data=True, 
                                      overwrite=True, 
                                      output_dir=path_tmp,
                                      file_name=None,
                                      )
    
            # unzip output
            gunzip(os.path.join(path_tmp,str(j)+"_def-img.nii.gz"))

        # merge final deformed time series
        data = nb.load(os.path.join(path_tmp,"0_def-img.nii"))
        data.header["dim"][4] = nt
        data_res = np.zeros(data.header["dim"][1:5])
        
        for j in range(nt):
            data_res[:,:,:,j] = nb.load(os.path.join(path_tmp,str(j)+"_def-img.nii")).get_fdata()            

        # time series path and basename
        path = os.path.dirname(input_epi[i])
        file = os.path.splitext(os.path.basename(input_epi[i]))[0]
        
        output = nb.Nifti1Image(data_res, data.affine, data.header)
        nb.save(output, os.path.join(path,"r"+file+"_linear.nii"))

        # delete intermediate files
        sh.rmtree(path_tmp, ignore_errors=True)

else:
    print("Number of time series and deformation are not the same!")
