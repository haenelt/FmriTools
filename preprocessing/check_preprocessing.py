"""
Spatial cross-correlation of the signal intensities in functional time series

This script computes the correlation between time steps in one or more concatenated functional time 
series. Only voxels within the brain are included in the analysis. The correlation strength between 
volumes is thought to be used as a qualitaty check for the realignment of time series.

created by Daniel Haenelt
Date created: 31-07-2019             
Last modified: 15-08-2019  
"""
import os
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.stats import pearsonr, shapiro

input = ["/home/daniel/Schreibtisch/processing/bsffp/run2/rdata.nii",
         "/home/daniel/Schreibtisch/processing/bsffp/run2/rdata.nii",
         ]
r_threshold = 0.8

""" do not edit below """

# font parameters for plots
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

# load data
for i in range(len(input)):
    if i == 0:
        data = nb.load(input[i]).get_fdata()
    else:
        data_temp = nb.load(input[i]).get_fdata()
        data = np.append(data, data_temp, axis=3)

# create output folder
path_output = os.path.join(os.path.dirname(input[0]),"correlation")
if not os.path.exists(path_output):
    os.makedirs(path_output)

# brain mask
mask = data[:,:,:,0].copy()
mask[mask < np.mean(mask)] = 0
mask[mask != 0] = 1

# time axis
t_shapiro = np.arange(0,np.shape(data)[3])
t_pearson = np.arange(2,np.shape(data)[3]+1)

# initialise arrays
r_pearson_0 = np.zeros_like(t_pearson, dtype="float")
r_pearson = np.zeros_like(t_pearson, dtype="float")
r_shapiro = np.zeros_like(t_shapiro, dtype="float")
p_pearson_0 = np.zeros_like(t_pearson, dtype="float")
p_pearson = np.zeros_like(t_pearson, dtype="float")
p_shapiro = np.zeros_like(t_shapiro, dtype="float")

# pool data from first volume
data_0 = data[:,:,:,0]
data_0[mask == 0] = np.nan
data_0 = data_0[~np.isnan(data_0)].flatten()

# shapiro-wilk of first volume
r_shapiro[0], p_shapiro[0] = shapiro(data_0)

# loop through time steps
for i in range(np.shape(data)[3]-1):

    data_i = data[:,:,:,i]
    data_j = data[:,:,:,i+1]
    
    data_i[mask == 0] = np.nan
    data_j[mask == 0] = np.nan
    
    data_i = data_i[~np.isnan(data_i)].flatten()
    data_j = data_j[~np.isnan(data_j)].flatten()

    r_pearson_0[i], p_pearson_0[i] = pearsonr(data_0,data_j)
    r_pearson[i], p_pearson[i] = pearsonr(data_i,data_j)
    r_shapiro[i+1], p_shapiro[i+1] = shapiro(data_j)

# percentage
res_pearson_0 = len(r_pearson_0[r_pearson_0 < r_threshold]) / len(data[0,0,0,:]) * 100
res_pearson = len(r_pearson[r_pearson < r_threshold]) / len(data[0,0,0,:]) * 100
res_shapiro = len(r_shapiro[r_shapiro < r_threshold]) / len(data[0,0,0,:]) * 100

# logfile
file = open(os.path.join(path_output,"correlation.txt"),"w")
file.write("Correlation threshold: "+str(r_threshold)+"\n")
file.write("pearson to volume 1: "+str(res_pearson_0)+"\n")
file.write("pearson to volume i-1: "+str(res_pearson)+"\n")
file.write("shapiro of volume i: "+str(res_shapiro)+"\n")
file.write("----------\n")
file.write("percentage of volumes below threshold")
file.close()

# plots
fig, ax = plt.subplots()
ax.plot(t_shapiro,r_shapiro, "r")
ax.set_xlabel("Time in TR")
ax.set_ylabel("r-value (Shapiro-Wilk test)")
ax.set_title("Shapiro-Wilk test for each time step")
ax.hlines(r_threshold,t_shapiro[0],t_shapiro[-1],linestyle="dashed")
fig.savefig(os.path.join(path_output,"r_shapiro"), format='svg', bbox_inches='tight')
plt.show()

fig, ax = plt.subplots()
ax.plot(t_shapiro,p_shapiro, "r")
ax.set_xlabel("Time in TR")
ax.set_ylabel("p-value (Shapiro-Wilk test)")
ax.set_title("Shapiro-Wilk test for each time step")
fig.savefig(os.path.join(path_output,"p_shapiro"), format='svg', bbox_inches='tight')
plt.show()

fig, ax = plt.subplots()
ax.plot(t_pearson,r_pearson_0, "r")
ax.set_xlabel("Time in TR")
ax.set_ylabel("r-value (Pearson correlation)")
ax.set_title("Time series spatial correlation to volume 1")
ax.hlines(r_threshold,t_pearson[0],t_pearson[-1],linestyle="dashed")
fig.savefig(os.path.join(path_output,"r_pearson_0"), format='svg', bbox_inches='tight')
plt.show()

fig, ax = plt.subplots()
ax.plot(t_pearson,p_pearson_0, "r")
ax.set_xlabel("Time in TR")
ax.set_ylabel("p-value (Pearson correlation)")
ax.set_title("Time series spatial correlation to volume 1")
fig.savefig(os.path.join(path_output,"p_pearson_0"), format='svg', bbox_inches='tight')
plt.show()

fig, ax = plt.subplots()
ax.plot(t_pearson,r_pearson, "r")
ax.set_xlabel("Time in TR")
ax.set_ylabel("r-value (Pearson correlation)")
ax.set_title("Time series spatial correlation to volume i-1")
ax.hlines(r_threshold,t_pearson[0],t_pearson[-1],linestyle="dashed")
fig.savefig(os.path.join(path_output,"r_pearson"), format='svg', bbox_inches='tight')
plt.show()

fig, ax = plt.subplots()
ax.plot(t_pearson,p_pearson, "r")
ax.set_xlabel("Time in TR")
ax.set_ylabel("p-value (Pearson correlation)")
ax.set_title("Time series spatial correlation to volume i-1")
fig.savefig(os.path.join(path_output,"p_pearson"), format='svg', bbox_inches='tight')
plt.show()