"""
Percent signal change

This scripts calculates the percent signal change between two conditions of a block design for a 
session consisting of several runs. First, a baseline correction of each time series is applied if 
not done before (i.e., if no file with prefix b is found). From the condition file which has to be 
in the SPM compatible *.mat format, time points for both blocks are defined. Time series for the 
whole time series (baseline) and all conditions can be converted to z-score. The percent signal 
change is computed as the difference of the mean between both conditions divided by the time series 
mean. The mean percent signal change of the whole session is taken as the average across single 
runs. The percent signal change is computed for both contrasts. If the outlier input array is not
empty, outlier volumes are discarded from the analysis.

created by Daniel Haenelt
Date created: 06-12-2018             
Last modified: 16-09-2019  
"""
import os
import datetime
import numpy as np
import nibabel as nb
from scipy.stats import zscore
from lib.processing import get_onset_vols

# input data
img_input = ["/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_1/udata.nii",
             "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_2/udata.nii",
             "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_3/udata.nii",
             "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_4/udata.nii",
             "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_5/udata.nii",
             "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_6/udata.nii",
             "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_7/udata.nii",
             "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_8/udata.nii",
             "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_9/udata.nii",
             "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_10/udata.nii",
             "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_11/udata.nii",
             "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_12/udata.nii",
             ]

cond_input = ["/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_1/logfiles/p8_SE_EPI1_Run1_odc_Cond.mat",
              "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_2/logfiles/p8_SE_EPI1_Run2_odc_Cond.mat",
              "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_3/logfiles/p8_SE_EPI1_Run3_odc_Cond.mat",
              "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_4/logfiles/p8_SE_EPI1_Run4_odc_Cond.mat",
              "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_5/logfiles/p8_SE_EPI1_Run5_odc_Cond.mat",
              "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_6/logfiles/p8_SE_EPI1_Run6_odc_Cond.mat",
              "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_7/logfiles/p8_SE_EPI1_Run7_odc_Cond.mat",
              "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_8/logfiles/p8_SE_EPI1_Run8_odc_Cond.mat",
              "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_9/logfiles/p8_SE_EPI1_Run9_odc_Cond.mat",
              "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_10/logfiles/p8_SE_EPI1_Run10_odc_Cond.mat",
              "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_11/logfiles/p8_SE_EPI1_Run11_odc_Cond.mat",
              "/data/pt_01880/V2STRIPES/p8/odc/SE_EPI1/Run_12/logfiles/p8_SE_EPI1_Run12_odc_Cond.mat",
              ]

outlier_input = []

# path to SPM12 folder
pathSPM = "/data/pt_01880/source/spm12"
pathLIB1 = "/home/raid2/haenelt/projects/scripts/lib/preprocessing"
pathLIB2 = "/home/raid2/haenelt/projects/scripts/lib/processing"

# parameters
TR = 3 # repetition time in s
cutoff_highpass = 180 # cutoff in s for baseline correction
skipvol = 3 # skip number of volumes in each block
condition1 = "left"
condition2 = "right"
name_output = "GE_EPI1"
use_z_score = False
use_lowpass = False
cutoff_lowpass = 0
order_lowpass = 0

""" do not edit below """

# prepare path and filename
path = []
file = []
for i in range(len(img_input)):
    path.append(os.path.split(img_input[i])[0])
    file.append(os.path.split(img_input[i])[1])

# output folder is taken from the first entry of the input list
path_output = os.path.join(os.path.dirname(os.path.dirname(path[0])),"results","percent","native")
if not os.path.exists(path_output):
    os.makedirs(path_output)

# get image header information from first entry of the input list
data_img = nb.load(img_input[0])
data_img.header["dim"][0] = 3
data_img.header["dim"][4] = 1
header = data_img.header
affine = data_img.affine

# get image dimension
dim = data_img.header["dim"][1:4]

mean_percent_signal1 = np.zeros(dim)
mean_percent_signal2 = np.zeros(dim)
for i in range(len(path)):
    
    if len(outlier_input) > 0:
        onsets1, onsets2 = get_onset_vols(cond_input[i], outlier_input[i], condition1, condition2, TR, skipvol)
    else:
        onsets1, onsets2 = get_onset_vols(cond_input[i], outlier_input, condition1, condition2, TR, skipvol)

    # lopass filter time series
    if use_lowpass:
        os.chdir(pathLIB2)
        os.system("matlab" + \
                  " -nodisplay -nodesktop -r " + \
                  "\"lowpass_filter(\'{0}\', {1}, {2}, {3}, \'{4}\'); exit;\"". \
                  format(img_input[i], TR, cutoff_lowpass, order_lowpass, pathSPM))
        
        # change input to lowpass filtered time series
        img_input[i] = os.path.join(path[i],"l"+file[i])

    # look for baseline corrected time series
    if not os.path.isfile(os.path.join(path[i],"b"+file[i])):
        os.chdir(pathLIB1)
        os.system("matlab" + \
                  " -nodisplay -nodesktop -r " + \
                  "\"baseline_correction(\'{0}\', {1}, {2}, \'{3}\'); exit;\"". \
                  format(img_input[i], TR, cutoff_highpass, pathSPM))

    # open baseline corrected data
    data_img = nb.load(os.path.join(path[i],"b"+file[i]))
    data_array = data_img.get_fdata()
    
    # sort volumes to conditions
    data_condition1 = data_array[:,:,:,onsets1]
    data_condition2 = data_array[:,:,:,onsets2]
    
    # z-score
    if use_z_score:
        data_condition1 = zscore(data_condition1, axis=3)
        data_condition2 = zscore(data_condition2, axis=3)
        data_array = zscore(data_array, axis=3)
    
    # mean
    data_condition1_mean = np.mean(data_condition1, axis=3)
    data_condition2_mean = np.mean(data_condition2, axis=3)
    data_baseline_mean = np.mean(data_array, axis=3)
    data_baseline_mean[data_baseline_mean == 0] = np.nan
    
    # percent signal change
    percent_signal1 = ( data_condition1_mean - data_condition2_mean ) / data_baseline_mean * 100
    percent_signal2 = ( data_condition2_mean - data_condition1_mean ) / data_baseline_mean * 100
    percent_signal1[np.isnan(percent_signal1)] = 0
    percent_signal2[np.isnan(percent_signal2)] = 0

    # sum volumes for each run
    mean_percent_signal1 += percent_signal1
    mean_percent_signal2 += percent_signal2
    
# divide by number of runs
mean_percent_signal1 = mean_percent_signal1 / len(path)
mean_percent_signal2 = mean_percent_signal2 / len(path)
    
# write output
output = nb.Nifti1Image(mean_percent_signal1, affine, header)
fileOUT = os.path.join(path_output,"percent_"+name_output+"_"+condition1+"_"+condition2+".nii")
nb.save(output,fileOUT)

output = nb.Nifti1Image(mean_percent_signal2, affine, header)
fileOUT = os.path.join(path_output,"percent_"+name_output+"_"+condition2+"_"+condition1+".nii")
nb.save(output,fileOUT)

# write log
fileID = open(os.path.join(path_output,"percent_info.txt"),"a")
fileID.write("script executed: "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
fileID.write("basename: "+name_output+"\n")
fileID.write("condition1: "+condition1+"\n")
fileID.write("condition2: "+condition2+"\n")
fileID.write("TR: "+str(TR)+"\n")
fileID.write("cutoff_highpass: "+str(cutoff_highpass)+"\n")
fileID.write("skipvol: "+str(skipvol)+"\n")
fileID.close()