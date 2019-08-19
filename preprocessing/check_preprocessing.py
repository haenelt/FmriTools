"""
Spatial cross-correlation of the signal intensities in functional time series

This script computes the correlation between time steps in one or more concatenated functional time 
series. Only voxels within the brain are included in the analysis. The correlation strength between 
volumes is thought to be used as a qualitaty check for the realignment of time series.

created by Daniel Haenelt
Date created: 31-07-2019             
Last modified: 19-08-2019  
"""
import os
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.stats import pearsonr, shapiro

input = [
        "/data/pt_01880/Experiment2_Rivalry/p2/odc/GE_EPI1/Run_1/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p2/odc/GE_EPI1/Run_2/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p2/odc/GE_EPI1/Run_3/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p2/odc/GE_EPI1/Run_4/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p2/odc/GE_EPI1/Run_5/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p2/odc/GE_EPI1/Run_6/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p2/odc/GE_EPI1/Run_7/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p2/odc/GE_EPI1/Run_8/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p2/odc/GE_EPI1/Run_9/udata.nii",
        "/data/pt_01880/Experiment2_Rivalry/p2/odc/GE_EPI1/Run_10/udata.nii",
        ]
r_threshold = 0.9

""" do not edit below """

# font parameters for plots
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

# create output folder
if len(input) < 2:
    path_output = os.path.join(os.path.dirname(input[0]),"correlation")
else:
    path_output = os.path.join(os.path.dirname(os.path.dirname(input[0])),"correlation")

if not os.path.exists(path_output):
    os.makedirs(path_output)

# pool data from first volume
data_0 = nb.load(input[0]).get_fdata()[:,:,:,0]

# open logfile
file = open(os.path.join(path_output,"correlation.txt"),"w")
file.write("Percentage of volumes below threshold\n\n")
file.write("Correlation threshold: "+str(r_threshold)+"\n")

r_pearson_0 = []
r_pearson = []
r_shapiro = []
p_pearson_0 = []
p_pearson = []
p_shapiro = []
neighbourFlag = 0 # omit correlation between neighbours for last time step
for i in range(len(input)):
    
    # load time series
    data_temp = nb.load(input[i]).get_fdata()   
    
    # print progress
    print("Time series "+str(i)+"/"+str(len(input)))
    
    npearson_0 = 0
    nshapiro = 0
    for j in range(np.shape(data_temp)[3]):
        
        # load first and current time step
        data_temp_0 = data_0.copy()
        data_temp_j = data_temp[:,:,:,j].copy()
        
        # get temporary mask
        mask_0 = data_temp_0.copy()
        mask_0[mask_0 < np.mean(mask_0[mask_0 > 0])] = 0
        mask_0[mask_0 != 0] = 1
        mask_j = data_temp_j.copy()
        mask_j[mask_j < np.mean(mask_j[mask_j > 0])] = 0
        mask_j[mask_j != 0] = 1        
        mask_0j = mask_0 * mask_j
        
        # mask time step
        data_temp_0 = data_temp_0[mask_0j == 1].flatten()
        data_temp_j = data_temp_j[mask_0j == 1].flatten()        
        
        # shapiro wilk
        r_tmp, p_tmp = shapiro(data_temp_j)
        r_shapiro = np.append(r_shapiro, r_tmp)
        p_shapiro = np.append(p_shapiro, p_tmp)

        if r_tmp < r_threshold:
            nshapiro += 1 
        
        # pearson correlation to first time step
        r_tmp, p_tmp = pearsonr(data_temp_0, data_temp_j)
        r_pearson_0 = np.append(r_pearson_0, r_tmp)
        p_pearson_0 = np.append(p_pearson_0, p_tmp)
        
        if r_tmp < r_threshold:
            npearson_0 += 1
        
        if j < np.shape(data_temp)[3]-1 and i < len(input):
            data_temp_1 = data_temp[:,:,:,j]
            data_temp_2 = data_temp[:,:,:,j+1]
        elif j == np.shape(data_temp)[3]-1 and i < len(input)-1:
            data_temp_1 = data_temp[:,:,:,j]
            data_temp_2 = nb.load(input[i+1]).get_fdata()[:,:,:,0]
        else:
            neighbourFlag = 1

        if neighbourFlag == 0:
            
            mask_1 = data_temp_1.copy()
            mask_1[mask_1 < np.mean(mask_1[mask_1 > 0])] = 0
            mask_1[mask_1 != 0] = 1
            mask_2 = data_temp_2.copy()
            mask_2[mask_2 < np.mean(mask_2[mask_2 > 0])] = 0
            mask_2[mask_2 != 0] = 1
            mask_12 = mask_1 * mask_2
    
            data_temp_1 = data_temp_1[mask_12 == 1].flatten()
            data_temp_2 = data_temp_2[mask_12 == 1].flatten()
        
            r_tmp, p_tmp = pearsonr(data_temp_1, data_temp_2)
            r_pearson = np.append(r_pearson, r_tmp)
            p_pearson = np.append(p_pearson, p_tmp)
               
    # percentage
    res_pearson_0 = npearson_0 / np.shape(data_temp)[3] * 100
    res_shapiro = nshapiro / np.shape(data_temp)[3] * 100
    
    # update logfile
    file.write("Run: "+str(i)+"\n")
    file.write("----------\n")
    file.write("pearson to volume 1: "+str(res_pearson_0)+"\n")
    file.write("shapiro of volume i: "+str(res_shapiro)+"\n\n\n")

# close logfile
file.close()

# save variables
np.savez(os.path.join(path_output,"correlation_data"),
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
fig.savefig(os.path.join(path_output,"r_shapiro"), format='svg', bbox_inches='tight')
#plt.show()

fig, ax = plt.subplots()
ax.plot(p_shapiro, "r")
ax.set_xlabel("Time in TR")
ax.set_ylabel("p-value (Shapiro-Wilk test)")
ax.set_title("Shapiro-Wilk test for each time step")
fig.savefig(os.path.join(path_output,"p_shapiro"), format='svg', bbox_inches='tight')
#plt.show()

fig, ax = plt.subplots()
ax.plot(r_pearson_0, "r")
ax.set_xlabel("Time in TR")
ax.set_ylabel("r-value (Pearson correlation)")
ax.set_title("Time series spatial correlation to volume 1")
ax.hlines(r_threshold,0,len(r_pearson_0),linestyle="dashed")
fig.savefig(os.path.join(path_output,"r_pearson_0"), format='svg', bbox_inches='tight')
#plt.show()

fig, ax = plt.subplots()
ax.plot(p_pearson_0, "r")
ax.set_xlabel("Time in TR")
ax.set_ylabel("p-value (Pearson correlation)")
ax.set_title("Time series spatial correlation to volume 1")
fig.savefig(os.path.join(path_output,"p_pearson_0"), format='svg', bbox_inches='tight')
#plt.show()

fig, ax = plt.subplots()
ax.plot(r_pearson, "r")
ax.set_xlabel("Time in TR")
ax.set_ylabel("r-value (Pearson correlation)")
ax.set_title("Time series spatial correlation to volume i-1")
ax.hlines(r_threshold,0,len(r_pearson),linestyle="dashed")
fig.savefig(os.path.join(path_output,"r_pearson"), format='svg', bbox_inches='tight')
#plt.show()

fig, ax = plt.subplots()
ax.plot(p_pearson, "r")
ax.set_xlabel("Time in TR")
ax.set_ylabel("p-value (Pearson correlation)")
ax.set_title("Time series spatial correlation to volume i-1")
fig.savefig(os.path.join(path_output,"p_pearson"), format='svg', bbox_inches='tight')
#plt.show()