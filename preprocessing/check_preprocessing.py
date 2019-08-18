"""
Spatial cross-correlation of the signal intensities in functional time series

This script computes the correlation between time steps in one or more concatenated functional time 
series. Only voxels within the brain are included in the analysis. The correlation strength between 
volumes is thought to be used as a qualitaty check for the realignment of time series.

created by Daniel Haenelt
Date created: 31-07-2019             
Last modified: 18-08-2019  
"""
import os
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.stats import pearsonr, shapiro

input = [
        "/data/pt_01880/preprocessing_test/method2/Run_1/rudata_linear_linear.nii",
        "/data/pt_01880/preprocessing_test/method2/Run_2/rudata_linear_linear.nii",
        ]
r_threshold = 0.95

""" do not edit below """

# font parameters for plots
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

# load first image in array
data = nb.load(input[0]).get_fdata()

# get brain mask from first volume
mask = data[:,:,:,0].copy()
mask[mask < np.mean(mask[mask > 0])] = 0
mask[mask != 0] = 1

# pool data from first volume
data_0 = data[:,:,:,0]
data_0[mask == 0] = np.nan
data_0 = data_0[~np.isnan(data_0)].flatten()

# create output folder
if len(input) < 2:
    path_output = os.path.join(os.path.dirname(input[0]),"correlation")
else:
    path_output = os.path.join(os.path.dirname(os.path.dirname(input[0])),"correlation")

if not os.path.exists(path_output):
    os.makedirs(path_output)

# initialise vectors
r_pearson_0 = []
r_pearson = []
r_shapiro = []
p_pearson_0 = []
p_pearson = []
p_shapiro = []

for i in range(len(input)):
    
    # load time series
    data_temp = nb.load(input[i]).get_fdata()
    
    # print iteration step
    print("Current time series: "+str(i))
    
    for j in range(np.shape(data_temp)[3]-1):
                
        # load time step
        data_temp_1 = data_temp[:,:,:,j]
        data_temp_2 = data_temp[:,:,:,j+1]
        data_temp_1[mask == 0] = np.nan
        data_temp_2[mask == 0] = np.nan
        data_temp_1 = data_temp_1[~np.isnan(data_temp_1)].flatten()
        data_temp_2 = data_temp_2[~np.isnan(data_temp_2)].flatten()
        
        # first volume of time series            
        if j == 0:
                
            r_tmp, p_tmp = shapiro(data_temp_1)
            r_shapiro = np.append(r_shapiro, r_tmp)
            p_shapiro = np.append(p_shapiro, p_tmp)
            
            r_tmp, p_tmp = pearsonr(data_0, data_temp_1)
            r_pearson_0 = np.append(r_pearson_0, r_tmp)
            p_pearson_0 = np.append(p_pearson_0, p_tmp)
            
        r_tmp, p_tmp = shapiro(data_temp_2)
        r_shapiro = np.append(r_shapiro, r_tmp)
        p_shapiro = np.append(p_shapiro, p_tmp)
        
        r_tmp, p_tmp = pearsonr(data_0, data_temp_2)
        r_pearson_0 = np.append(r_pearson_0, r_tmp)
        p_pearson_0 = np.append(p_pearson_0, p_tmp)
        
        r_tmp, p_tmp = pearsonr(data_temp_1, data_temp_2)
        r_pearson = np.append(r_pearson, r_tmp)
        p_pearson = np.append(p_pearson, p_tmp)

# percentage
res_pearson_0 = len(r_pearson_0[r_pearson_0 < r_threshold]) / len(r_pearson_0) * 100
res_pearson = len(r_pearson[r_pearson < r_threshold]) / len(r_pearson) * 100
res_shapiro = len(r_shapiro[r_shapiro < r_threshold]) / len(r_shapiro) * 100

# logfile
file = open(os.path.join(path_output,"correlation.txt"),"w")
file.write("Correlation threshold: "+str(r_threshold)+"\n")
file.write("pearson to volume 1: "+str(res_pearson_0)+"\n")
file.write("pearson to volume i-1: "+str(res_pearson)+"\n")
file.write("shapiro of volume i: "+str(res_shapiro)+"\n")
file.write("----------\n")
file.write("percentage of volumes below threshold")
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
