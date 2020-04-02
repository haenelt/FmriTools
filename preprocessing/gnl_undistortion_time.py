"""
Gradient nonlinearity correction of time series

This scripts calls the HCP toolbox to correct for gradient nonlinearities in the input time series.

Before running the script, login to queen via ssh and set the fsl environment by calling FSL in 
the terminal. The corrected time series gets the suffix _gnlcorr.

created by Daniel Haenelt
Date created: 27-01-2020
Last modified: 28-01-2020
"""
import os
import shutil as sh
import numpy as np
import nibabel as nb
from nipype.interfaces.fsl import ExtractROI
from nipype.interfaces.fsl.preprocess import ApplyWarp
from lib.io.get_filename import get_filename
from lib.preprocessing.gnl_correction import gnl_correction

# input
file_in = [
    "/data/pt_01880/Experiment3_Stripes/p3/colour/GE_EPI1/Run_1/data.nii",
    ]

file_bash = "/data/hu_haenelt/projects/gradunwarp/apply_grad.sh"
file_coeff = "/data/hu_haenelt/projects/gradunwarp/7t_coeff.grad"
python3_env = "daniel"
python2_env = "daniel2"
cleanup = True

""" do not edit below """

for i in range(len(file_in)):
    
    # get fileparts of input
    path_file, name_file, ext_file = get_filename(file_in[i])
    
    # filenames
    file_vol0 = os.path.join(path_file,name_file+"_vol0"+ext_file)
    file_out = os.path.join(path_file,name_file+"_gnlcorr"+ext_file)
    
    # extract first volume
    fslroi = ExtractROI()
    fslroi.inputs.in_file = file_in[i]
    fslroi.inputs.roi_file = file_vol0
    fslroi.inputs.output_type = "NIFTI"
    fslroi.inputs.t_min = 0
    fslroi.inputs.t_size = 1
    fslroi.run()

    # exexute gnl correction
    gnl_correction(file_vol0, 
                   file_bash, 
                   file_coeff, 
                   python3_env, 
                   python2_env, 
                   path_file, 
                   False)
    
    # apply warp to first volume
    applywarp = ApplyWarp()
    applywarp.inputs.in_file = file_in[i]
    applywarp.inputs.ref_file = file_in[i]
    applywarp.inputs.relwarp = True
    applywarp.inputs.field_file = os.path.join(path_file,"grad","warp.nii.gz")
    applywarp.inputs.output_type = "NIFTI"
    applywarp.inputs.out_file = file_out
    applywarp.inputs.interp = "spline"
    applywarp.run()

    # normalise warped output image to initial intensity range
    data_img = nb.load(file_in[i])
    data_array = data_img.get_fdata()
    max_data = np.max(data_array)
    min_data = np.min(data_array)
    
    data_img = nb.load(file_out)
    data_array = data_img.get_fdata()
    data_array[data_array < min_data] = 0
    data_array[data_array > max_data] = max_data
    
    output = nb.Nifti1Image(data_array, data_img.affine, data_img.header)
    nb.save(output, file_out)
    
    # remove vol0
    os.remove(file_vol0)
    os.remove(os.path.join(path_file,name_file+"_vol0_gnlcorr"+ext_file))
    
    # clean intermediate files
    if cleanup:
        sh.rmtree(os.path.join(path_file,"grad"), ignore_errors=True)
