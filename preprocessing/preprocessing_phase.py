"""
Contrast-enhanced mean epi

This scripts preprocesses the phase data from an epi time series to get the phase variance over
time and an mean epi magnitude image with enhanced GM/WM contrast. The following steps are executed:
    1. motion correction to magnitude time series
    2. unwrap phase time series
    3. apply correction to phase
    4. get mean magnitude and phase
    5. get phase variance
    6. rescale mean phase image
    7. filter mean phase image
    8. weight mean magnitude epi

Before running the script, login to queen via ssh and set the afni environment by calling AFNI in 
the terminal.

created by Daniel Haenelt
Date created: 10-01-2020
Last modified: 20-04-2020
"""
import os
import numpy as np
import nibabel as nb
from nipype.interfaces import afni
from nighres.intensity import phase_unwrapping
from lib.io import get_filename
from lib.utils import get_mean, get_std
from lib.preprocessing import deweight_mask

# input data
input_magn = "/data/pt_01880/test_data/data/data.nii"
input_phase = "/data/pt_01880/test_data/data/data_phase.nii"

# parameters
phase_max = 0.25 # threshold mean phase data
std_max = 0.25 # threshold for mask generation before phase filtering
sigma_gaussian = 10.0 # sigma for gaussian filter


""" do not edit below """

path_magn, name_magn, ext_magn = get_filename(input_magn)
path_phase, name_phase, ext_phase = get_filename(input_phase)

# make subfolders
path_moco = os.path.join(path_magn,'moco')
if not os.path.exists(path_moco):
    os.makedirs(path_moco)
    
"""
motion correction of magnitude data
"""    
moco = afni.Volreg()
moco.inputs.in_file = os.path.join(path_magn,name_magn+ext_magn)
moco.inputs.out_file = os.path.join(path_magn,'u'+name_magn+ext_magn)
moco.inputs.args = '-base 0 -twopass -float -clipit'
moco.inputs.zpad = 4
moco.inputs.interp = 'Fourier'
moco.inputs.outputtype = 'NIFTI'
moco.inputs.md1d_file = os.path.join(path_moco,'max_disp.1D')
moco.inputs.oned_file = os.path.join(path_moco,'moco_params.1D')
moco.inputs.oned_matrix_save = os.path.join(path_moco,'moco_matrix.1D')
moco.run()

"""
unwrap phase data
"""
magn_img = nb.load(input_magn)
phase_img = nb.load(input_phase)

# change header information
magn_img.header["datatype"] = 64
phase_img.header["dim"][0] = 3
phase_img.header["dim"][4] = 1

# unwrap single time steps independently
phase_array = phase_img.get_fdata()
phase_array2 = np.zeros_like(phase_array)
for i in range(np.shape(phase_array)[3]):
    print("Unwrap volume "+str(i))
    temp_in = nb.Nifti1Image(phase_array[:,:,:,i], phase_img.affine, phase_img.header)
    temp_out = phase_unwrapping(temp_in, mask=None, nquadrants=3, tv_flattening=True, tv_scale=0.5,
                                save_data=False, overwrite=False, output_dir=False, file_name=False)
    
    phase_array2[:,:,:,i] = temp_out["result"].get_fdata()

# save unwrapped time series
name_phase = name_phase+"_unwrap"
output = nb.Nifti1Image(phase_array2, magn_img.affine, magn_img.header)
nb.save(output, os.path.join(path_phase, name_phase+ext_phase))

"""
apply motion correction to unwrapped phase time series
"""
allineate = afni.Allineate()
allineate.inputs.in_file = os.path.join(path_phase, name_phase+ext_phase)
allineate.inputs.out_file = os.path.join(path_phase, "u"+name_phase+ext_phase)
allineate.inputs.in_matrix= os.path.join(path_moco,'moco_matrix.1D')
allineate.inputs.outputtype = 'NIFTI'
allineate.inputs.final_interpolation = 'linear'
allineate.inputs.no_pad = True # do not use zero-padding on the base image
allineate.run()

"""
get mean time series
"""
name_magn = "u"+name_magn
name_phase = "u"+name_phase
get_mean(os.path.join(path_magn, name_magn+ext_magn), path_magn, name_magn, type="mean")
get_mean(os.path.join(path_phase, name_phase+ext_phase), path_phase, name_phase, type="mean")

"""
phase variance over time
"""
get_std(os.path.join(path_phase, name_phase+ext_phase), path_phase, name_phase)

"""
rescale phase data
"""
phase_img = nb.load(os.path.join(path_phase,"mean_"+name_phase+ext_phase))
phase_array = phase_img.get_fdata()

# threshold phase data
phase_array[phase_array > phase_max*np.max(phase_array)] = phase_max*np.max(phase_array) 
phase_array[phase_array < -phase_max*np.max(phase_array)] = -phase_max*np.max(phase_array) 

# rescale
phase_range = np.max(phase_array) - np.min(phase_array)
phase_array = (phase_array - np.min(phase_array)) / phase_range

# invert contrast
phase_array = 1 - phase_array

output = nb.Nifti1Image(phase_array, phase_img.affine, phase_img.header)
nb.save(output, os.path.join(path_phase,"mean_"+name_phase+ext_phase))

"""
enhanced contrast
"""
magn_img = nb.load(os.path.join(path_magn,"mean_"+name_magn+ext_magn))
magn_array = magn_img.get_fdata()

# filter phase
phase_array = deweight_mask(os.path.join(path_phase,"mean_"+name_phase+ext_phase), 
                            os.path.join(path_phase,"std_"+name_phase+ext_phase), 
                            std_max, 
                            sigma_gaussian, 
                            write_output=None, 
                            path_output=None)

# weight magnitude image
magn_array = magn_array * phase_array

output = nb.Nifti1Image(magn_array, magn_img.affine, magn_img.header)
nb.save(output, os.path.join(path_magn,"mean_"+name_magn+"_enhanced"+ext_magn))
