"""
Contrast-to-noise ratio of functional time series (task-based activation)

This scripts calculates the contrast-to-noise ratio (CNR) in percent from functional time series  
containing task-based activation following a block design. The input can be a list of several runs. 
From the condition file which has to be in the SPM compatible *.mat format, time points for both an
experimental and a baseline condition are extracted. CNR is computed as absolute difference between 
both conditions divided by the standard deviation of the baseline condition. The CNR of the whole 
session is taken as the average across single runs. Similar computations of CNR can be found in 
Scheffler et al. (2016). If the outlier input array is not empty, outlier volumes are discarded from 
the analysis. Optionally, the time series can be filtered by a highpass filter. The input images 
should be in nifti format.

Before running the script, login to queen via ssh and set the afni environment by calling AFNI in 
the terminal.

created by Daniel Haenelt
Date created: 03-05-2019             
Last modified: 17-03-2020  
"""
import os
import datetime
import numpy as np
import nibabel as nb
from lib.io.get_filename import get_filename
from lib.processing import get_onset_vols

# input data
img_input = [
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_1/budata.nii",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_2/budata.nii",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_3/budata.nii",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_4/budata.nii",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_5/budata.nii",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_6/budata.nii",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_7/budata.nii",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_8/budata.nii",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_9/budata.nii",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_10/budata.nii",
        ]

cond_input = [
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_1/logfiles/p3_GE_EPI2_Run1_odc_Cond.mat",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_2/logfiles/p3_GE_EPI2_Run2_odc_Cond.mat",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_3/logfiles/p3_GE_EPI2_Run3_odc_Cond.mat",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_4/logfiles/p3_GE_EPI2_Run4_odc_Cond.mat",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_5/logfiles/p3_GE_EPI2_Run5_odc_Cond.mat",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_6/logfiles/p3_GE_EPI2_Run6_odc_Cond.mat",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_7/logfiles/p3_GE_EPI2_Run7_odc_Cond.mat",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_8/logfiles/p3_GE_EPI2_Run8_odc_Cond.mat",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_9/logfiles/p3_GE_EPI2_Run9_odc_Cond.mat",
        "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_10/logfiles/p3_GE_EPI2_Run10_odc_Cond.mat",
        ]

outlier_input = []

# parameters
condition0 = "rest" # baseline condition
condition1 = "left" # experimental condition
TR = 3 # repetition time in s
skip_vol = 2 # skip number of volumes in each block
use_highpass = False
cutoff_highpass = 180 # cutoff in s for baseline correction
name_sess = "GE_EPI2"
name_output = "super"

# path to SPM12 folder
pathSPM = "/data/pt_01880/source/spm12"
pathLIB = "/home/raid2/haenelt/projects/scripts/lib/preprocessing"

""" do not edit below """

# get path from first entry
path_file, _, _ = get_filename(img_input[0])

# make output folder
path_output = os.path.join(os.path.dirname(os.path.dirname(path_file)),"results","cnr","native")
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

mean_cnr = np.zeros(dim)
for i in range(len(img_input)):

    # get filename
    path_file, name_file, ext_file = get_filename(img_input[i])
    
    # get condition specific onsets
    onsets0 = get_onset_vols(cond_input[i], outlier_input[i], condition0, TR, skip_vol)
    onsets1 = get_onset_vols(cond_input[i], outlier_input[i], condition1, TR, skip_vol)
    
    # highpass filter time series
    if use_highpass:
        os.chdir(pathLIB)
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
    data_condition0 = data_array[:,:,:,onsets0]
    data_condition1 = data_array[:,:,:,onsets1]
       
    # mean
    data_condition0_mean = np.mean(data_condition0, axis=3)
    data_condition1_mean = np.mean(data_condition1, axis=3)
    data_condition0_std = np.std(data_condition0, axis=3)
    data_condition0_std[data_condition0_std == 0] = np.nan
    
    # percent signal change
    cnr = ( np.abs(data_condition1_mean - data_condition0_mean) ) / data_condition0_std * 100
    cnr[np.isnan(cnr)] = 0

    # sum volumes for each run
    mean_cnr += cnr
    
# divide by number of runs
mean_cnr /= len(img_input)

# name of output files
if len(name_output) and len(name_sess):
    fileOUT = os.path.join(path_output,"cnr_"+name_output+"_"+condition1+"_"+condition0+"_"+name_sess+".nii")
elif len(name_output) and not len(name_sess):
    fileOUT = os.path.join(path_output,"cnr_"+name_output+"_"+condition1+"_"+condition0+".nii")
elif not len(name_output) and len(name_sess):
    fileOUT = os.path.join(path_output,"cnr_"+condition1+"_"+condition0+"_"+name_sess+".nii")
else:
    fileOUT = os.path.join(path_output,"cnr_"+condition1+"_"+condition0+".nii")

# write output
output = nb.Nifti1Image(mean_cnr, affine, header)
nb.save(output,fileOUT)

# write log
fileID = open(os.path.join(path_output,"cnr_info.txt"),"a")
fileID.write("script executed: "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
fileID.write("session: "+name_sess+"\n")
fileID.write("basename: "+name_output+"\n")
fileID.write("condition0: "+condition0+"\n")
fileID.write("condition1: "+condition1+"\n")
fileID.write("TR: "+str(TR)+"\n")
fileID.write("skip_vol: "+str(skip_vol)+"\n")
fileID.write("highpass: "+str(use_highpass)+"\n")
fileID.write("cutoff highpass: "+str(cutoff_highpass)+"\n")
fileID.close()