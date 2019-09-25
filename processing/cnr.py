"""
Contrast-to-noise ratio of functional time series (task-based activation)

This scripts calculates the contrast-to-noise ratio (CNR) in percent from functional time series  
containing task-based activation following a block design. The input can be a list of several runs. 
First, a baseline correction of each time series is applied if not done before (i.e., if no file 
with prefix b is found). From the condition file which has to be in the SPM compatible *.mat format, 
time points for both conditions are defined. CNR is computed as absolute difference between both 
conditions divided by the standard deviation of the second condition. The second condition should be 
a baseline condition if the CNR should make any sense. The CNR of the whole session is taken as the 
average across single runs. Similar computations of CNR can be found in Scheffler et al. (2016). If
the outlier input array is not empty, outlier volumes are discarded from the analysis. Optionally 
(if n != 0), the time series can be upsampled. The input images should be in nifti format. 

Before running the script, login to queen via ssh and set the afni environment by calling AFNI in 
the terminal.

created by Daniel Haenelt
Date created: 03-05-2019             
Last modified: 25-09-2019  
"""
import os
import datetime
import numpy as np
import nibabel as nb
from lib.processing import get_onset_vols
from lib.utils import upsample_time_series

# input data
img_input = [
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_1/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_2/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_3/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_4/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_5/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_6/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_7/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_8/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_9/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_10/udata.nii",
        ]

cond_input = [
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_1/logfiles/p3_GE_EPI2_Run1_nonrivalry_Cond.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_2/logfiles/p3_GE_EPI2_Run2_nonrivalry_Cond.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_3/logfiles/p3_GE_EPI2_Run3_nonrivalry_Cond.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_4/logfiles/p3_GE_EPI2_Run4_nonrivalry_Cond.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_5/logfiles/p3_GE_EPI2_Run5_nonrivalry_Cond.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_6/logfiles/p3_GE_EPI2_Run6_nonrivalry_Cond.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_7/logfiles/p3_GE_EPI2_Run7_nonrivalry_Cond.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_8/logfiles/p3_GE_EPI2_Run8_nonrivalry_Cond.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_9/logfiles/p3_GE_EPI2_Run9_nonrivalry_Cond.mat",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_10/logfiles/p3_GE_EPI2_Run10_nonrivalry_Cond.mat",
        ]

outlier_input = [
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_1/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_2/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_3/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_4/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_5/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_6/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_7/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_8/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_9/logfiles/outlier_regressor.txt",
        "/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI2/Run_10/logfiles/outlier_regressor.txt",
        ]

# path to SPM12 folder
pathSPM = "/data/pt_01880/source/spm12"
pathLIB = "/home/raid2/haenelt/projects/scripts/lib/preprocessing"

# parameters
TR = 3 # repetition time in s
cutoff_highpass = 180 # cutoff in s for baseline correction
skipvol = 0 # skip number of volumes in each block
condition1 = "right"
condition2 = "rest"
name_sess = "GE_EPI2"
name_output = ""
n = 0 # upsampling factor

""" do not edit below """

# change to lib folder
os.chdir(pathLIB)

# prepare path and filename
path = []
file = []
for i in range(len(img_input)):
    path.append(os.path.split(img_input[i])[0])
    file.append(os.path.splitext(os.path.split(img_input[i])[1])[0])

# output folder is taken from the first entry of the input list
path_output = os.path.join(os.path.dirname(os.path.dirname(path[0])),"results","cnr","native")
if not os.path.exists(path_output):
    os.makedirs(path_output)

# get TR for upsampled time series
if n:
    TR = TR / n

# get upsampled time series
if n:
    for i in range(len(img_input)):
        file[i] = file[i]+"_upsampled"
        if not os.path.isfile(os.path.join(path[i],file[i]+".nii")):
            upsample_time_series(img_input[i], n)

# get image header information from first entry of the input list
data_img = nb.load(os.path.join(path[0],file[0]+".nii"))
data_img.header["dim"][0] = 3
data_img.header["dim"][4] = 1
header = data_img.header
affine = data_img.affine

# get image dimension
dim = data_img.header["dim"][1:4]

mean_cnr = np.zeros(dim)
for i in range(len(path)):
    
    if len(outlier_input) > 0:
        onsets1, onsets2 = get_onset_vols(cond_input[i], outlier_input[i], condition1, condition2, TR, skipvol)
    else:
        onsets1, onsets2 = get_onset_vols(cond_input[i], outlier_input, condition1, condition2, TR, skipvol)

    # look for baseline corrected time series
    if not os.path.isfile(os.path.join(path[i],"b"+file[i]+".nii")):
        os.system("matlab" + \
                  " -nodisplay -nodeskop -r " + \
                  "\"baseline_correction(\'{0}\', {1}, {2}, \'{3}\'); exit;\"". \
                  format(os.path.join(path[i],file[i]+".nii"), TR, cutoff_highpass, pathSPM))

    # open baseline corrected data
    data_img = nb.load(os.path.join(path[i],"b"+file[i]+".nii"))
    data_array = data_img.get_fdata()
    
    # sort volumes to conditions
    data_condition1 = data_array[:,:,:,onsets1]
    data_condition2 = data_array[:,:,:,onsets2]
       
    # mean
    data_condition1_mean = np.mean(data_condition1, axis=3)
    data_condition2_mean = np.mean(data_condition2, axis=3)
    data_condition2_std = np.std(data_condition2, axis=3)
    data_condition2_std[data_condition2_std == 0] = np.nan
    
    # percent signal change
    cnr = ( np.abs(data_condition1_mean - data_condition2_mean) ) / data_condition2_std * 100
    cnr[np.isnan(cnr)] = 0

    # sum volumes for each run
    mean_cnr += cnr
    
# divide by number of runs
mean_cnr = mean_cnr / len(path)

# name of output files
if len(name_output) and len(name_sess):
    fileOUT = os.path.join(path_output,"cnr_"+name_output+"_"+condition1+"_"+condition2+"_"+name_sess+".nii")
elif len(name_output) and not len(name_sess):
    fileOUT = os.path.join(path_output,"cnr_"+name_output+"_"+condition1+"_"+condition2+".nii")
elif not len(name_output) and len(name_sess):
    fileOUT = os.path.join(path_output,"cnr_"+condition1+"_"+condition2+"_"+name_sess+".nii")
else:
    fileOUT = os.path.join(path_output,"cnr_"+condition1+"_"+condition2+".nii")

# write output
output = nb.Nifti1Image(mean_cnr, affine, header)
nb.save(output,fileOUT)

# write log
fileID = open(os.path.join(path_output,"cnr_info.txt"),"a")
fileID.write("script executed: "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
fileID.write("session: "+name_sess+"\n")
fileID.write("basename: "+name_output+"\n")
fileID.write("condition1: "+condition1+"\n")
fileID.write("condition2: "+condition2+"\n")
fileID.write("TR: "+str(TR)+"\n")
fileID.write("cutoff_highpass: "+str(cutoff_highpass)+"\n")
fileID.write("skipvol: "+str(skipvol)+"\n")
fileID.write("n: "+str(n)+"\n")
fileID.close()