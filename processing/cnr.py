"""
Percent signal change

This scripts calculates the percent signal change between two conditions of a block design for a 
session consisting of several runs. First, a baseline correction of each time series is applied if 
not done before (i.e., if no file with prefix b is found). From the condition file which has to be 
in the SPM compatible *.mat format, time points for both blocks are defined. Time series for the 
whole time series (baseline) and all conditions can be converted to z-score. The percent signal 
change is computed as the difference of the mean between both conditions divided by the time series 
mean. The mean percent signal change of the whole session is taken as the average across single 
runs. The percent signal change is computed for both contrasts.

created by Daniel Haenelt
Date created: 06-12-2018             
Last modified: 09-01-2019  
"""
import sys
import os
import datetime
import numpy as np
import nibabel as nb
from scipy.io import loadmat
from scipy.stats import zscore

# input data
img_input = ["/nobackup/actinium2/haenelt/VasoTest/flicker/Run_1/uvaso_basis_corrected.nii",
             "/nobackup/actinium2/haenelt/VasoTest/flicker/Run_2/uvaso_basis_corrected.nii",
             "/nobackup/actinium2/haenelt/VasoTest/flicker/Run_3/uvaso_basis_corrected.nii",
             "/nobackup/actinium2/haenelt/VasoTest/flicker/Run_4/uvaso_basis_corrected.nii",
             "/nobackup/actinium2/haenelt/VasoTest/flicker/Run_5/uvaso_basis_corrected.nii",
             "/nobackup/actinium2/haenelt/VasoTest/flicker/Run_6/uvaso_basis_corrected.nii",
             ]

cond_input = ["/nobackup/actinium2/haenelt/VasoTest/flicker/Run_1/logfiles/VasoTest_flicker_Run1_Cond.mat",
              "/nobackup/actinium2/haenelt/VasoTest/flicker/Run_2/logfiles/VasoTest_flicker_Run2_Cond.mat",
              "/nobackup/actinium2/haenelt/VasoTest/flicker/Run_3/logfiles/VasoTest_flicker_Run3_Cond.mat",
              "/nobackup/actinium2/haenelt/VasoTest/flicker/Run_4/logfiles/VasoTest_flicker_Run4_Cond.mat",
              "/nobackup/actinium2/haenelt/VasoTest/flicker/Run_5/logfiles/VasoTest_flicker_Run5_Cond.mat",
              "/nobackup/actinium2/haenelt/VasoTest/flicker/Run_6/logfiles/VasoTest_flicker_Run6_Cond.mat",
              ]

# path to SPM12 folder
pathSPM = "/data/pt_01880/source/spm12"
pathLIB = "/home/raid2/haenelt/projects/scripts/lib/preprocessing"

# parameters
TR = 2.5 # repetition time in s
cutoff_highpass = 100 # cutoff in s for baseline correction
skipvol = 4 # skip number of volumes in each block
condition1 = "on"
condition2 = "off"
name_output = "vaso"
use_z_score = False

""" do not edit below """

# change to lib folder
os.chdir(pathLIB)

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

# output filenames
name_sess = os.path.basename(os.path.dirname(path[0]))

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
    
    # load condition file
    cond = loadmat(cond_input[i])

    # get condition information
    names = np.concatenate(np.concatenate(cond["names"]))
    onsets = np.concatenate(cond["onsets"])
    durations = np.concatenate(np.concatenate(np.concatenate(cond["durations"])))

    # check if condition names exist
    if not condition1 in names:
        sys.exit("condition1 not found in the condition_file")

    if not condition2 in names:
        sys.exit("condition2 not found in the condition_file")

    # index of both conditions    
    c1 = np.where(condition1 == names)[0][0]
    c2 = np.where(condition2 == names)[0][0]

    # onsets for both conditions
    onsets1 = np.round(onsets[c1][0] / TR + skipvol).astype(int)
    onsets2 = np.round(onsets[c2][0] / TR + skipvol).astype(int)

    # durations for both conditions
    durations1 = np.round(durations[c1] / TR - skipvol).astype(int)
    durations2 = np.round(durations[c2] / TR - skipvol).astype(int)

    # sort all volumes to be considered in both conditions
    temp = onsets1.copy()
    for j in range(durations1 - 1):
        onsets1 = np.append(onsets1, temp + j + 1)
    onsets1 = np.sort(onsets1)

    temp = onsets2.copy()
    for j in range(durations2 - 1):
        onsets2 = np.append(onsets2, temp + j + 1)
    onsets2 = np.sort(onsets2)

    # look for baseline corrected time series
    if not os.path.isfile(os.path.join(path[i],"b"+file[i])):
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
    percent_signal1[np.isnan(percent_signal1)] = 0

    # sum volumes for each run
    mean_percent_signal1 += percent_signal1
    mean_percent_signal2 += percent_signal2
    
# divide by number of runs
mean_percent_signal1 = mean_percent_signal1 / len(path)
mean_percent_signal2 = mean_percent_signal2 / len(path)
    
# write output
output = nb.Nifti1Image(mean_percent_signal1, affine, header)
fileOUT = os.path.join(path_output,"percent_"+name_output+"_"+condition1+"_"+condition2+"_"+name_sess+".nii")
nb.save(output,fileOUT)

output = nb.Nifti1Image(mean_percent_signal2, affine, header)
fileOUT = os.path.join(path_output,"percent_"+name_output+"_"+condition2+"_"+condition1+"_"+name_sess+".nii")
nb.save(output,fileOUT)

# write log
fileID = open(os.path.join(path_output,"percent_info.txt"),"a")
fileID.write("script executed: "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
fileID.write("basename: "+name_output+"\n")
fileID.write("data set: "+name_sess+"\n")
fileID.write("condition1: "+condition1+"\n")
fileID.write("condition2: "+condition2+"\n")
fileID.write("TR: "+str(TR)+"\n")
fileID.write("cutoff_highpass: "+str(cutoff_highpass)+"\n")
fileID.write("skipvol: "+str(skipvol)+"\n")
fileID.close()