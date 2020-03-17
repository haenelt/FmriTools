"""
Percent signal change

This scripts calculates the percent signal change between two conditions of a block design for a 
session consisting of several runs. First, a baseline correction of each time series is applied if 
not done before (i.e., if no file with prefix b is found). From the condition file which has to be 
in the SPM compatible *.mat format, time points for both experimental blocks and baseline blocks are 
defined. The percent signal change (psc) is computed as the difference of the mean between both
conditions relative to either the mean of the whole time series or to a baseline condition. The mean
psc of the whole session is taken as the average across single runs. PSC is computed for both
contrast directions. If the outlier input array is not epty, outlier volumes are discarded from the
analysis. The input images should be in nifti format.

Before running the script, login to queen via ssh and set the afni environment by calling AFNI in 
the terminal.

created by Daniel Haenelt
Date created: 06-12-2018         
Last modified: 17-03-2020
"""
import os
import datetime
import numpy as np
import nibabel as nb
from scipy.stats import zscore
from lib.io.get_filename import get_filename
from lib.processing import get_onset_vols

# input data
img_input = [
        "/data/pt_01880/analysis/VASO3/mean_budata.nii",
        ]

cond_input = [
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_1/logfiles/p3_VASO3_Run1_odc_Cond.mat",
        ]

outlier_input = []

# path to SPM12 folder
pathSPM = "/data/pt_01880/source/spm12"
pathLIB1 = "/home/raid2/haenelt/projects/scripts/lib/preprocessing"
pathLIB2 = "/home/raid2/haenelt/projects/scripts/lib/processing"

# parameters
TR = 2.5 # repetition time in s
cutoff_highpass = 180 # cutoff in s for baseline correction
skip_vol = 2 # skip number of volumes in each block
condition0 = "rest"
condition1 = "right"
condition2 = "rest"
name_sess = "VASO3"
name_output = ""
use_z_score = False
use_lowpass = False
baseline_all = False
cutoff_lowpass = 0
order_lowpass = 0

""" do not edit below """

# prepare path and filename
path = []
file = []
for i in range(len(img_input)):
    path_temp, file_temp, _ = get_filename(img_input[i])
    path.append(path_temp)
    file.append(file_temp)

# output folder is taken from the first entry of the input list
path_output = os.path.join(os.path.dirname(os.path.dirname(path[0])),"results","raw","native")
if not os.path.exists(path_output):
    os.makedirs(path_output)

# outlier
if not len(outlier_input):
    outlier_input = np.zeros(len(img_input))
    

# get image header information from first entry of the input list
data_img = nb.load(os.path.join(path[0],file[0]+".nii"))
data_img.header["dim"][0] = 3
data_img.header["dim"][4] = 1
header = data_img.header
affine = data_img.affine

# get image dimension
dim = data_img.header["dim"][1:4]

mean_percent_signal1 = np.zeros(dim)
mean_percent_signal2 = np.zeros(dim)
for i in range(len(path)):
    
    # get condition specific onsets
    onsets0 = get_onset_vols(cond_input[i], outlier_input[i], condition0, TR, skip_vol)
    onsets1 = get_onset_vols(cond_input[i], outlier_input[i], condition1, TR, skip_vol)
    onsets2 = get_onset_vols(cond_input[i], outlier_input[i], condition2, TR, skip_vol)

    # lopass filter time series
    if use_lowpass:
        os.chdir(pathLIB2)
        os.system("matlab" + \
                  " -nodisplay -nodesktop -r " + \
                  "\"lowpass_filter(\'{0}\', {1}, {2}, {3}, \'{4}\'); exit;\"". \
                  format(os.path.join(path[i],file[i]+".nii"), TR, cutoff_lowpass, order_lowpass, pathSPM))
        
        # change input to lowpass filtered time series
        file[i] = "l"+file[i]

    # look for baseline corrected time series
    if not os.path.isfile(os.path.join(path[i],"b"+file[i]+".nii")):
        os.chdir(pathLIB1)
        os.system("matlab" + \
                  " -nodisplay -nodesktop -r " + \
                  "\"baseline_correction(\'{0}\', {1}, {2}, \'{3}\'); exit;\"". \
                  format(os.path.join(path[i],file[i]+".nii"), TR, cutoff_highpass, pathSPM))

    # open baseline corrected data
    data_img = nb.load(os.path.join(path[i],"b"+file[i]+".nii"))
    data_array = data_img.get_fdata()
    
    # sort volumes to conditions
    data_condition0 = data_array[:,:,:,onsets0]
    data_condition1 = data_array[:,:,:,onsets1]
    data_condition2 = data_array[:,:,:,onsets2]
    
    # z-score
    if use_z_score:
        data_condition0 = zscore(data_condition0, axis=3)
        data_condition1 = zscore(data_condition1, axis=3)
        data_condition2 = zscore(data_condition2, axis=3)
        data_array = zscore(data_array, axis=3)
    
    # mean
    data_condition0_mean = np.mean(data_condition0, axis=3)
    data_condition1_mean = np.mean(data_condition1, axis=3)
    data_condition2_mean = np.mean(data_condition2, axis=3)
    data_baseline_mean = np.mean(data_array, axis=3)
    data_baseline_mean[data_baseline_mean == 0] = np.nan
    
    # percent signal change
    if baseline_all:
        percent_signal1 = ( data_condition1_mean - data_condition2_mean ) / data_baseline_mean * 100
        percent_signal2 = ( data_condition2_mean - data_condition1_mean ) / data_baseline_mean * 100
    else:
        percent_signal1 = ( data_condition1_mean - data_condition2_mean ) / data_condition0_mean * 100
        percent_signal2 = ( data_condition2_mean - data_condition1_mean ) / data_condition0_mean * 100

    percent_signal1[np.isnan(percent_signal1)] = 0
    percent_signal2[np.isnan(percent_signal2)] = 0

    # sum volumes for each run
    mean_percent_signal1 += percent_signal1
    mean_percent_signal2 += percent_signal2
    
# divide by number of runs
mean_percent_signal1 = mean_percent_signal1 / len(path)
mean_percent_signal2 = mean_percent_signal2 / len(path)

# name of output files
if len(name_output) and len(name_sess):
    fileOUT1 = os.path.join(path_output,"psc_"+name_output+"_"+condition1+"_"+condition2+"_"+name_sess+".nii")
    fileOUT2 = os.path.join(path_output,"psc_"+name_output+"_"+condition2+"_"+condition1+"_"+name_sess+".nii")
elif len(name_output) and not len(name_sess):
    fileOUT1 = os.path.join(path_output,"psc_"+name_output+"_"+condition1+"_"+condition2+".nii")
    fileOUT2 = os.path.join(path_output,"psc_"+name_output+"_"+condition2+"_"+condition1+".nii")
elif not len(name_output) and len(name_sess):
    fileOUT1 = os.path.join(path_output,"psc_"+condition1+"_"+condition2+"_"+name_sess+".nii")
    fileOUT2 = os.path.join(path_output,"psc_"+condition2+"_"+condition1+"_"+name_sess+".nii")
else:
    fileOUT1 = os.path.join(path_output,"psc_"+condition1+"_"+condition2+".nii")
    fileOUT2 = os.path.join(path_output,"psc_"+condition2+"_"+condition1+".nii")

# write output
output = nb.Nifti1Image(mean_percent_signal1, affine, header)
nb.save(output,fileOUT1)

output = nb.Nifti1Image(mean_percent_signal2, affine, header)
nb.save(output,fileOUT2)

# write log
fileID = open(os.path.join(path_output,"psc_info.txt"),"a")
fileID.write("script executed: "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
fileID.write("session: "+name_sess+"\n")
fileID.write("basename: "+name_output+"\n")
fileID.write("condition1: "+condition1+"\n")
fileID.write("condition2: "+condition2+"\n")
fileID.write("TR: "+str(TR)+"\n")
fileID.write("cutoff_highpass: "+str(cutoff_highpass)+"\n")
fileID.write("skip_vol: "+str(skip_vol)+"\n")
fileID.close()