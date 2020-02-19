"""
Spatial cross-correlation of the signal intensities in functional time series

This scripts intends to identify corrupted volumes and runs in a set of functional time series to 
get a quantifiable and reproducible exclusion criterion (cf. Marquardt et al. 2017; Bergmann et al.
2019). Reference volumes can be for example the session mean. A regressor of no interest is written 
for each input time series to denote volumes below threshold.

created by Daniel Haenelt
Date created: 31-07-2019             
Last modified: 19-02-2020  
"""
import os
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, shapiro
from lib.io.get_filename import get_filename

input = [
    "/data/pt_01880/temp/run1/udata.nii",
    "/data/pt_01880/temp/run2/udata.nii",
    ]
input_ref = "/data/pt_01880/temp/input/mean_epi.nii"
input_mask_ref = "/data/pt_01880/temp/input/mask_def-img3.nii.gz"
r_threshold = 0.95

""" do not edit below """

# get filename from first input entry
_, name_file, _ = get_filename(input[0])

# create output folder
if len(input) < 2:
    path_output = os.path.join(os.path.dirname(input[0]),"correlation")
else:
    path_output = os.path.join(os.path.dirname(os.path.dirname(input[0])),"correlation")

if not os.path.exists(path_output):
    os.makedirs(path_output)

# load reference volume
data_0 = nb.load(input_ref).get_fdata()

if len(input_mask_ref) > 0:
    mask_0 = nb.load(input_mask_ref).get_fdata()

# open logfile
file = open(os.path.join(path_output,"correlation_"+name_file+".txt"),"w")
file.write("Percentage of volumes below threshold\n")
file.write("Correlation threshold: "+str(r_threshold)+"\n\n")

r_pearson_0 = []
r_pearson = []
r_shapiro = []
p_pearson_0 = []
p_pearson = []
p_shapiro = []
for i in range(len(input)):
    
    # load time series
    data_temp = nb.load(input[i]).get_fdata()   
    
    # print progress
    print("Time series "+str(i+1)+"/"+str(len(input)))
    
    # create output folder for regressor of no interest
    path_logfile = os.path.join(os.path.dirname(input[i]),"logfiles")
    
    if not os.path.exists(path_logfile):
        os.makedirs(path_logfile)
    
    # open logfile for regressor of no interest
    file2 = open(os.path.join(path_logfile,"correlation_regressor_"+name_file+".txt"),"w")
    
    pearson_run = 0
    npearson_0 = 0
    nshapiro = 0
    for j in range(np.shape(data_temp)[3]):
        
        # load reference and current time step
        data_temp_0 = data_0.copy()
        data_temp_j = data_temp[:,:,:,j].copy()
        
        # mask time step
        if len(input_mask_ref) > 0:
            data_temp_0 = data_temp_0[mask_0 == 1].flatten()
            data_temp_j = data_temp_j[mask_0 == 1].flatten()       
        else:
            data_temp_0 = data_temp_0.flatten()
            data_temp_j = data_temp_j.flatten()       
        
        # shapiro wilk
        r_tmp, p_tmp = shapiro(data_temp_j)
        r_shapiro = np.append(r_shapiro, r_tmp)
        p_shapiro = np.append(p_shapiro, p_tmp)

        if r_tmp < r_threshold:
            nshapiro += 1 
        
        # pearson correlation to reference
        r_tmp, p_tmp = pearsonr(data_temp_0, data_temp_j)
        r_pearson_0 = np.append(r_pearson_0, r_tmp)
        p_pearson_0 = np.append(p_pearson_0, p_tmp)
        
        # sum for within-run correlation
        pearson_run += r_tmp
        
        if r_tmp < r_threshold:
            npearson_0 += 1
            file2.write("1\n")
        else:
            file2.write("0\n")
        
        if j < np.shape(data_temp)[3]-1 and i < len(input):
            data_temp_1 = data_temp[:,:,:,j].copy()
            data_temp_2 = data_temp[:,:,:,j+1].copy()
            
        elif j == np.shape(data_temp)[3]-1 and i < len(input)-1:
            data_temp_1 = data_temp[:,:,:,j].copy()
            data_temp_2 = nb.load(input[i+1]).get_fdata()[:,:,:,0]
            
        else:
            break
            
        # mask time step
        if len(input_mask_ref) > 0:
            data_temp_1 = data_temp_1[mask_0 == 1].flatten()
            data_temp_2 = data_temp_2[mask_0 == 1].flatten() 
        else:
            data_temp_1 = data_temp_1.flatten()
            data_temp_2 = data_temp_2.flatten() 
        
        r_tmp, p_tmp = pearsonr(data_temp_1, data_temp_2)
        r_pearson = np.append(r_pearson, r_tmp)
        p_pearson = np.append(p_pearson, p_tmp)
               
    # close logfile for regressor of no interest
    file2.close()
    
    # percentage below threshold
    res_pearson_0 = npearson_0 / np.shape(data_temp)[3] * 100
    res_shapiro = nshapiro / np.shape(data_temp)[3] * 100
    
    # average within-run correlation
    pearson_run = pearson_run / np.shape(data_temp)[3]
    
    # update logfile
    file.write("Run: "+str(i+1)+"\n")
    file.write("----------\n")
    file.write("Pearson (average within run): "+str(pearson_run)+"\n")
    file.write("Outlier percentage (pearson to ref): "+str(res_pearson_0)+"\n")
    file.write("Outlier percentage (shapiro): "+str(res_shapiro)+"\n\n\n")

# close logfile
file.close()

# save variables
np.savez(os.path.join(path_output,"correlation_"+name_file),
         r_shapiro=r_shapiro, p_shapiro=p_shapiro,
         r_pearson=r_pearson, p_pearson=p_pearson,
         r_pearson_0=r_pearson_0, p_pearson_0=p_pearson_0,
         )

# plots
fig, ax = plt.subplots()
ax.plot(r_shapiro, "r")
ax.set_xlabel("Time in TR")
ax.set_ylabel("r-value (Shapiro-Wilk test)")
ax.set_title("Shapiro-Wilk test for each time step")
ax.hlines(r_threshold,0,len(r_shapiro),linestyle="dashed")
fig.savefig(os.path.join(path_output,"r_shapiro_"+name_file+".png"), format='png', bbox_inches='tight')
#plt.show()

fig, ax = plt.subplots()
ax.plot(p_shapiro, "r")
ax.set_xlabel("Time in TR")
ax.set_ylabel("p-value (Shapiro-Wilk test)")
ax.set_title("Shapiro-Wilk test for each time step")
fig.savefig(os.path.join(path_output,"p_shapiro_"+name_file+".png"), format='png', bbox_inches='tight')
#plt.show()

fig, ax = plt.subplots()
ax.plot(r_pearson_0, "r")
ax.set_xlabel("Time in TR")
ax.set_ylabel("r-value (Pearson correlation)")
ax.set_title("Time series spatial correlation to volume ref")
ax.hlines(r_threshold,0,len(r_pearson_0),linestyle="dashed")
fig.savefig(os.path.join(path_output,"r_pearson_0_"+name_file+".png"), format='png', bbox_inches='tight')
#plt.show()

fig, ax = plt.subplots()
ax.plot(p_pearson_0, "r")
ax.set_xlabel("Time in TR")
ax.set_ylabel("p-value (Pearson correlation)")
ax.set_title("Time series spatial correlation to volume ref")
fig.savefig(os.path.join(path_output,"p_pearson_0_"+name_file+".png"), format='png', bbox_inches='tight')
#plt.show()

fig, ax = plt.subplots()
ax.plot(r_pearson, "r")
ax.set_xlabel("Time in TR")
ax.set_ylabel("r-value (Pearson correlation)")
ax.set_title("Time series spatial correlation to volume i-1")
ax.hlines(r_threshold,0,len(r_pearson),linestyle="dashed")
fig.savefig(os.path.join(path_output,"r_pearson_"+name_file+".png"), format='png', bbox_inches='tight')
#plt.show()

fig, ax = plt.subplots()
ax.plot(p_pearson, "r")
ax.set_xlabel("Time in TR")
ax.set_ylabel("p-value (Pearson correlation)")
ax.set_title("Time series spatial correlation to volume i-1")
fig.savefig(os.path.join(path_output,"p_pearson_"+name_file+".png"), format='png', bbox_inches='tight')
#plt.show()