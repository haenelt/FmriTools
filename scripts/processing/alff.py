# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import nibabel as nb

# local inputs
from fmri_tools.preprocessing.get_nuisance_mask import get_nuisance_mask
from fmri_tools.preprocessing.get_nuisance_regressor import get_nuisance_regressor
from fmri_tools.processing.get_alff import get_alff


"""
ALFF

The purpose of the following script is to compute ALFF and fALFF for a 
resting-state time series. The time series is baseline corrected and nuisance 
regressors are taken either from mean WM and CSF masks or based on biopac 
recordings. Optionally, no nuisance regression is performed.

The script needs an installation of freesurfer, ants adn afni.

created by Daniel Haenelt
Date created: 02-03-2019
Last modified: 13-10-2020
"""

# input
anatomy = "" # T1w full brain anatomy (e.g. orig)
function = "/data/pt_01983/func/resting_state2/uadata.nii" # baseline uncorrected
deformation = "" # deformation ana -> epi
biopac_input = "" # *.mat file
path_output = "/data/pt_01983/func/resting_state2/alff/native"

# parameters
TR = 3 # repetition time in s
cutoff_highpass = 270 # cutoff frequency for baseline correction in 1/Hz
nerode_wm = 1 # number of wm mask eroding iterations
nerode_csf = 1 # number of csf mask eroding iterations
hp_freq = 0.01 # highpass cutoff frequency (bandpass filter) in Hz
lp_freq = 0.08 # lowpass cutoff frequency (bandpass filter) in Hz

# analysis type
nuisance_regression = False
segmentation = False
biopac = False
cleanup = True

# do not edit below

# make output folder
if not os.path.exists(path_output):
    os.makedirs(path_output)

# get path and filenames
path = os.path.dirname(function)
file = os.path.basename(function)
bfile = "b" + file # filename of baseline corrected time series
rfile = "r" + file # filename of residual time series

# physiological noise regression
if nuisance_regression:
    
    # baseline correction
    os.system("matlab" + \
              " -nodisplay -nodesktop -r " + \
              "\"baseline_correction(\'{0}\', {1}, {2}); exit;\"". \
              format(function, TR, cutoff_highpass))

    if biopac:
        
        # get biopac regressors
        os.system("matlab" + \
                  " -nodisplay -nodesktop -r " + \
                  "\"get_biopac_regressor(\'{0}\', \'{1}\', \'{2}\', {3}); exit;\"". \
                  format(os.path.join(path,bfile), biopac_input, path_output, TR))
    
    else:

        # get wm and csf mask
        get_nuisance_mask(anatomy, deformation, path_output, nerode_wm, 
                          nerode_csf, segmentation, cleanup)
    
        # set mask to zero where function is equal to zero
        func_array = nb.load(function).get_fdata()
        func_array = func_array[:,:,:,0]
    
        wm = nb.load(os.path.join(path_output,"wm_mask.nii.gz"))
        wm_array = wm.get_fdata()
        wm_array[func_array == 0] = 0
        output = nb.Nifti1Image(wm_array, wm.affine, wm.header)
        nb.save(output,os.path.join(path_output,"wm_mask.nii.gz"))

        csf = nb.load(os.path.join(path_output,"csf_mask.nii.gz"))
        csf_array = csf.get_fdata()
        csf_array[func_array == 0] = 0
        output = nb.Nifti1Image(csf_array, csf.affine, csf.header)
        nb.save(output,os.path.join(path_output,"csf_mask.nii.gz"))

        # get nuisance regressor
        get_nuisance_regressor(os.path.join(path,bfile), 
                               os.path.join(path_output,"wm_mask.nii.gz"), 
                               os.path.join(path_output,"csf_mask.nii.gz"), 
                               path_output)
    
    # nuisance regression
    if cleanup:
        clean_glm = 1
    else:
        clean_glm = 0

    os.system("matlab" + \
              " -nodisplay -nodesktop -r " + \
              "\"regress_physio(\'{0}\', \'{1}\', {2}, {3}, \'{4}\', {5}); exit;\"". \
              format(function, 
                     os.path.join(path_output,"nuisance_regressor.txt"), 
                     TR, cutoff_highpass, path_output, clean_glm))

    # get alff    
    get_alff(os.path.join(path_output,rfile), TR, path_output, hp_freq, lp_freq, cleanup)
    
    # remove baseline corrected time series
    if cleanup:
        os.remove(os.path.join(path,bfile))

else:

    # get alff
    get_alff(function, TR, path_output, hp_freq, lp_freq, cleanup)
