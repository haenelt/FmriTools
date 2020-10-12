# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb
from nipype.interfaces import afni
from nighres.intensity import phase_unwrapping

# local inputs
from fmri_tools.io import get_filename
from fmri_tools.utils import get_mean, get_std
from fmri_tools.preprocessing import deweight_mask


"""
Contrast-enhanced mean epi

This scripts preprocesses the phase data from an epi time series to get the 
phase variance over time and an mean epi magnitude image with enhanced GM/WM 
contrast. The following steps are executed:
    1. motion correction to magnitude time series
    2. unwrap phase time series
    3. apply correction to phase
    4. remove outlier volumes from timeseries
    5. get mean magnitude and phase
    6. get phase variance
    7. rescale mean phase image
    8. filter mean phase image
    9. weight mean magnitude epi

Before running the script, login to queen via ssh and set the afni environment 
by calling AFNI in the terminal.

created by Daniel Haenelt
Date created: 10-01-2020
Last modified: 12-10-2020
"""

# input data
input_magn = "/data/pt_01880/temp_phase/data.nii"
input_phase = "/data/pt_01880/temp_phase/data_phase.nii"

# parameters
phase_max = 0.25 # threshold mean phase data
std_max = 0.25 # threshold for mask generation before phase filtering
sigma_gaussian = 10.0 # sigma for gaussian filter
outlier_params = [0.4, 0.8, 0.5, 1.0, 2.0] # mm short, mm long, deg short, deg long, z

# add spm and fmri_tools to path
pathSPM = "/data/pt_01880/source/spm12"
pathFMRITOOLS = "/data/hu_haenelt/projects/FmriTools/fmri_tools"

# do not edit below

# change to preprocessing folder in fmri_tools
os.chdir(os.path.join(pathFMRITOOLS,"preprocessing"))

path_magn, name_magn, ext_magn = get_filename(input_magn)
path_phase, name_phase, ext_phase = get_filename(input_phase)

# make subfolders
path_moco = os.path.join(path_magn,'diagnosis')
if not os.path.exists(path_moco):
    os.makedirs(path_moco)
    
# motion correction of magnitude data
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

# plot motion parameters
os.system("matlab" + \
          " -nodisplay -nodesktop -r " + \
          "\"plot_moco(\'{0}\', \'afni\', \'{1}\', \'{2}\'); exit;\"". \
          format(os.path.join(path_moco,'moco_params.1D'), path_moco, 'rp'))

# sum motion outliers
os.system("matlab" + \
          " -nodisplay -nodesktop -r " + \
          "\"get_outlier(\'{0}\', \'{1}\', \'{2}\', \'afni\', \'{3}\', \'{4}\'); exit;\"". \
          format(os.path.join(path_moco,"moco_params.1D"), 
                 os.path.join(path_magn,'u'+name_magn+ext_magn), 
                 outlier_params, 
                 os.path.join(path_magn,"outlier"),
                 pathSPM))
    
# unwrap phase data
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

# apply motion correction to unwrapped phase time series
allineate = afni.Allineate()
allineate.inputs.in_file = os.path.join(path_phase, name_phase+ext_phase)
allineate.inputs.out_file = os.path.join(path_phase, "u"+name_phase+ext_phase)
allineate.inputs.in_matrix= os.path.join(path_moco,'moco_matrix.1D')
allineate.inputs.outputtype = 'NIFTI'
allineate.inputs.final_interpolation = 'linear'
allineate.inputs.no_pad = True # do not use zero-padding on the base image
allineate.run()

# set new basenames to target realigned and temporary timeseries
name_magn = "u"+name_magn
name_phase = "u"+name_phase
name_magn_temp = name_magn+"_temp"
name_phase_temp = name_phase+"_temp"

# remove outliers
t = np.loadtxt(os.path.join(path_magn,"outlier", "outlier_regressor_"+name_magn+".txt")).astype(int)
    
# remove vols from magnitude data
magn_img = nb.load(os.path.join(path_magn, name_magn+ext_magn))
magn_array = magn_img.get_data()
magn_array = magn_array[:,:,:,t==0]
output = nb.Nifti1Image(magn_array, magn_img.affine, magn_img.header)
nb.save(output,os.path.join(path_magn, name_magn_temp+ext_magn))
    
# remove vols from phase data
phase_img = nb.load(os.path.join(path_phase, name_phase+ext_phase))
phase_array = phase_img.get_data()
phase_array = phase_array[:,:,:,t==0]
output = nb.Nifti1Image(phase_array, phase_img.affine, phase_img.header)
nb.save(output,os.path.join(path_phase, name_phase_temp+ext_phase))

# get mean time series
get_mean(os.path.join(path_magn, name_magn_temp+ext_magn), path_magn, name_magn, type="mean")
get_mean(os.path.join(path_phase, name_phase_temp+ext_phase), path_phase, name_phase, type="mean")

# phase variance over time
get_std(os.path.join(path_phase, name_phase_temp+ext_phase), path_phase, name_phase)

# rescale phase data
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

# enhanced contrast
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

# remove copied timeseries
os.remove(os.path.join(path_magn, name_magn_temp+ext_magn))
os.remove(os.path.join(path_phase, name_phase_temp+ext_phase))

# outlier percentage
outlier = np.loadtxt(os.path.join(path_magn,"outlier","outlier_regressor_"+name_magn+".txt"))
outlier_percentage = np.sum(outlier) / len(outlier) * 100

# write summary
fileID = open(os.path.join(path_moco,"preprocessing_summary_"+name_magn+".txt"),"w")
fileID.write("List of input parameters\n")
fileID.write("----------\n\n")
fileID.write("Preprocessed data\n")
fileID.write(input_magn+"\n")
fileID.write(input_phase+"\n")
fileID.write("threshold mean phase: "+str(phase_max)+"\n")
fileID.write("threshold std phase: "+str(std_max)+"\n")
fileID.write("gaussian filter (sigma): "+str(sigma_gaussian)+"\n\n")
fileID.write("Percentage of within-run motion and intensity outliers\n")
fileID.write("motion threshold (mm, short): "+str(outlier_params[0])+"\n")
fileID.write("motion threshold (mm, long): "+str(outlier_params[1])+"\n")
fileID.write("motion threshold (deg, short): "+str(outlier_params[2])+"\n")
fileID.write("motion threshold (deg, long): "+str(outlier_params[3])+"\n")
fileID.write("intensity threshold (z-score): "+str(outlier_params[4])+"\n")
fileID.write("----------\n\n")
fileID.write(str(round(outlier_percentage,2))+"\n")
fileID.close() 
