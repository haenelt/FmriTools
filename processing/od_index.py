"""
OD index

This scripts calculates a defined ocular dominance (OD) index for a session consisting of several
runs. From the condition file which has to be in the SPM compatible *.mat format, time points for 
both experimental conditions are extracted. Condition time points can be converted to z-scores. The 
OD index is computed by dividing each condition mean by the condition mean or max within a 
predefined mask before computing the difference of both conditions. The index for the whole session 
is taken as the average across single runs. If the outlier input array is not empty, outlier volumes 
are discarded from the analysis. Optionally, the time series can be filtered by a lowpass and a 
highpass filter. The input images should be in nifti format.

Before running the script, login to queen via ssh and set the afni environment by calling AFNI in 
the terminal.

created by Daniel Haenelt
Date created: 16-09-2019             
Last modified: 17-03-2020  
"""
import sys
import os
import datetime
import numpy as np
import nibabel as nb
from scipy.stats import zscore
from nighres.registration import apply_coordinate_mappings
from lib.io.get_filename import get_filename
from lib.processing import get_onset_vols

# input data
img_input = [
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_1/uadata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_2/uadata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_3/uadata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_4/uadata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_5/uadata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_6/uadata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_7/uadata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_8/uadata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_9/uadata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_10/uadata.nii",
        ]

cond_input = [
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_1/logfiles/p3_GE_EPI1_Run1_rivalry_Cond_threshold.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_2/logfiles/p3_GE_EPI1_Run2_rivalry_Cond_threshold.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_3/logfiles/p3_GE_EPI1_Run3_rivalry_Cond_threshold.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_4/logfiles/p3_GE_EPI1_Run4_rivalry_Cond_threshold.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_5/logfiles/p3_GE_EPI1_Run5_rivalry_Cond_threshold.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_6/logfiles/p3_GE_EPI1_Run6_rivalry_Cond_threshold.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_7/logfiles/p3_GE_EPI1_Run7_rivalry_Cond_threshold.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_8/logfiles/p3_GE_EPI1_Run8_rivalry_Cond_threshold.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_9/logfiles/p3_GE_EPI1_Run9_rivalry_Cond_threshold.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_10/logfiles/p3_GE_EPI1_Run10_rivalry_Cond_threshold.mat",
        ]

outlier_input = [
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_1/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_2/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_3/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_4/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_5/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_6/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_7/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_8/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_9/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI1/Run_10/logfiles/outlier_regressor.txt",
        ]

mask_input = "/data/pt_01880/Experiment2_Rivalry/p3/retinotopy/prf/vol/native/roi_target_sigma_threshold.nii"
epi2orig_input = "/data/pt_01880/Experiment2_Rivalry/p3/deformation/retinotopy/epi2orig.nii.gz"
orig2epi_input = "/data/pt_01880/Experiment2_Rivalry/p3/deformation/nonrivalry/GE_EPI2/orig2epi.nii.gz"

# parameters
condition1 = "left" # experimental condition 1
condition2 = "right" # experimental condition 2
TR = 3 # repetition time in s
skip_vol = 3 # skip number of volumes in each block
baseline_calculation = "mean" # mean or max
use_z_score = False
use_highpass = True
use_lowpass = False
cutoff_highpass = 180 # cutoff in s for baseline correction
cutoff_lowpass = 0
order_lowpass = 0
name_sess = "GE_EPI2"
name_output = ""

# path to SPM12 folder
pathSPM = "/data/pt_01880/source/spm12"
pathLIB1 = "/home/raid2/haenelt/projects/scripts/lib/preprocessing"
pathLIB2 = "/home/raid2/haenelt/projects/scripts/lib/processing"

""" do not edit below """

# get path from first entry
path_file, _, _ = get_filename(img_input[0])

# make output folder
path_output = os.path.join(os.path.dirname(os.path.dirname(path_file)),"results","od","native")
if not os.path.exists(path_output):
    os.makedirs(path_output)

# get image header information
data_img = nb.load(img_input[0])
data_img.header["dim"][0] = 3
data_img.header["dim"][4] = 1
header = data_img.header
affine = data_img.affine

# get image dimension
dim = data_img.header["dim"][1:4]

# get outlier dummy array if not outlier input
if not len(outlier_input):
    outlier_input = np.zeros(len(img_input))

# deform mask
apply_coordinate_mappings(mask_input, # input 
                          epi2orig_input, # cmap1
                          orig2epi_input,
                          interpolation = "nearest", # nearest or linear
                          padding = "zero", # closest, zero or max
                          save_data = True, # save output data to file (boolean)
                          overwrite = True, # overwrite existing results (boolean)
                          output_dir = path_output, # output directory
                          file_name = "mask" # base name with file extension for output
                          )

# rename final deformations
os.rename(os.path.join(path_output,"mask_def-img.nii.gz"),
          os.path.join(path_output,"mask.nii.gz"))

# load mask
mask = nb.load(os.path.join(path_output,"mask.nii.gz")).get_fdata()

mean_od_index1 = np.zeros(dim)
mean_od_index2 = np.zeros(dim)
for i in range(len(img_input)):
    
    # get filename
    path_file, name_file, ext_file = get_filename(img_input[i])
    
    # get condition specific onsets
    onsets1 = get_onset_vols(cond_input[i], outlier_input[i], condition1, TR, skip_vol)
    onsets2 = get_onset_vols(cond_input[i], outlier_input[i], condition2, TR, skip_vol)
    
    # lowpass filter time series
    if use_lowpass:
        os.chdir(pathLIB2)
        os.system("matlab" + \
                  " -nodisplay -nodesktop -r " + \
                  "\"lowpass_filter(\'{0}\', {1}, {2}, {3}, \'{4}\'); exit;\"". \
                  format(os.path.join(path_file,name_file+ext_file), TR, cutoff_lowpass, 
                         order_lowpass, pathSPM))
        
        # change input to lowpass filtered time series
        name_file = "l" + name_file

    # highpass filter time series
    if use_highpass:
        os.chdir(pathLIB1)
        os.system("matlab" + \
                  " -nodisplay -nodesktop -r " + \
                  "\"baseline_correction(\'{0}\', {1}, {2}, \'{3}\'); exit;\"". \
                  format(os.path.join(path_file,name_file+ext_file), TR, cutoff_highpass, pathSPM))

        # change input to highpass filtered time series
        name_file = "b" + name_file

    # open baseline corrected data
    data_img = nb.load(os.path.join(path_file,name_file+ext_file))
    data_array = data_img.get_fdata()
    
    # sort volumes to conditions
    data_condition1 = data_array[:,:,:,onsets1]
    data_condition2 = data_array[:,:,:,onsets2]
    
    # z-score
    if use_z_score:
        data_condition1 = zscore(data_condition1, axis=3)
        data_condition2 = zscore(data_condition2, axis=3)
    
    # mean
    data_condition1_mean = np.mean(data_condition1, axis=3)
    data_condition2_mean = np.mean(data_condition2, axis=3)
    
    if baseline_calculation is "mean":
        data_condition1_baseline = np.mean(data_condition1_mean[mask == 1])
        data_condition2_baseline = np.mean(data_condition2_mean[mask == 1])
    elif baseline_calculation is "max":
        data_condition1_baseline = np.max(data_condition1_mean[mask == 1])
        data_condition2_baseline = np.max(data_condition2_mean[mask == 1])   
    else:
        sys.exit("select either mean or max for baseline calculation!")
    
    # od index
    data_condition1_rel = data_condition1_mean / data_condition1_baseline
    data_condition2_rel = data_condition2_mean / data_condition2_baseline
        
    od_index1 = ( data_condition1_rel - data_condition2_rel ) * 100
    od_index2 = ( data_condition2_rel - data_condition1_rel ) * 100

    # sum volumes for each run
    mean_od_index1 += od_index1
    mean_od_index2 += od_index2
    
# divide by number of runs
mean_od_index1 /= len(img_input)
mean_od_index2 /= len(img_input)

# name of output files
if len(name_output) and len(name_sess):
    fileOUT1 = os.path.join(path_output,"od_"+name_output+"_"+condition1+"_"+condition2+"_"+name_sess+".nii")
    fileOUT2 = os.path.join(path_output,"od_"+name_output+"_"+condition2+"_"+condition1+"_"+name_sess+".nii")
elif len(name_output) and not len(name_sess):
    fileOUT1 = os.path.join(path_output,"od_"+name_output+"_"+condition1+"_"+condition2+".nii")
    fileOUT2 = os.path.join(path_output,"od_"+name_output+"_"+condition2+"_"+condition1+".nii")
elif not len(name_output) and len(name_sess):
    fileOUT1 = os.path.join(path_output,"od_"+condition1+"_"+condition2+"_"+name_sess+".nii")
    fileOUT2 = os.path.join(path_output,"od_"+condition2+"_"+condition1+"_"+name_sess+".nii")
else:
    fileOUT1 = os.path.join(path_output,"od_"+condition1+"_"+condition2+".nii")
    fileOUT2 = os.path.join(path_output,"od_"+condition2+"_"+condition1+".nii")

# write output
output = nb.Nifti1Image(mean_od_index1, affine, header)
nb.save(output,fileOUT1)

output = nb.Nifti1Image(mean_od_index2, affine, header)
nb.save(output,fileOUT2)

# write log
fileID = open(os.path.join(path_output,"od_info.txt"),"a")
fileID.write("script executed: "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
fileID.write("session: "+name_sess+"\n")
fileID.write("basename: "+name_output+"\n")
fileID.write("condition1: "+condition1+"\n")
fileID.write("condition2: "+condition2+"\n")
fileID.write("TR: "+str(TR)+"\n")
fileID.write("skip_vol: "+str(skip_vol)+"\n")
fileID.write("baseline calculation: "+baseline_calculation+"\n")
fileID.write("z-score: "+str(use_z_score)+"\n")
fileID.write("highpass: "+str(use_highpass)+"\n")
fileID.write("lowpass: "+str(use_lowpass)+"\n")
fileID.write("cutoff highpass: "+str(cutoff_highpass)+"\n")
fileID.write("cutoff lowpass: "+str(cutoff_lowpass)+"\n")
fileID.write("order lowpass: "+str(order_lowpass)+"\n")
fileID.close()