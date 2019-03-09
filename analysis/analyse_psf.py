"""
PSF estimation of multipolar phase-encoded paradigm.

TO DO:
    *fit error: covariance matrix
    *fit alternative: considering doing an exponential fit (Lorentzian) instead of a Gaussian

created by Daniel Haenelt
Date created: 08-03-2019
Last modified: 08-03-2019
"""
import os
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# input
path_ortho = "/data/pt_01880/V2STRIPES/p6/anatomy/ortho"
path_data = "/data/pt_01880/V2STRIPES/p6/psf/results/surf"
path_output = "/data/pt_01880/test"
file_patch = "occip4"

# parameters
multipol = [2, 4, 6, 8, 10, 12, 14] # conditions
layer = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9] # layer
min_distance = 10 # minimal cortical distance
max_distance = np.Inf # maximal cortical distance
min_percentile = 90 # minimum percentile
max_percentile = 100 # maximum percentile
layer_percentile = 0 # layer from which percentile is computed
bin_size = 20 # number of bins
hemi = ["lh", "rh"] # hemisphere

""" do not edit below """

# make output folder
if not os.path.exists(path_output):
    os.makedirs(path_output)

def get_percentile(input_data, input_cmap, input_distance, min_percentile, max_percentile, min_distance, max_distance):
    """
    This function selects the voxels which will be considered in further analysis. Voxels are 
    selected based on their corresponding minimum or maximum length and based on the chose 
    percentile of their coherence value. 
    """
    
    # load data
    data = nb.load(input_data).get_fdata()
    cmap = nb.load(input_cmap).get_fdata().astype(int)
    distance = nb.load(input_distance).get_fdata()
    
    # get coordinates within ROI
    distance = distance[cmap != 0]
    cmap = cmap[cmap != 0]
    cmap = cmap[np.all([~np.isnan(distance), distance != 0], axis=0)]
    distance = distance[np.all([~np.isnan(distance), distance != 0], axis=0)]
       
    # sort by distance
    cmap = cmap[np.all([distance > min_distance, distance < max_distance], axis=0)]
    distance = distance[np.all([distance > min_distance, distance < max_distance], axis=0)]    
    
    # sort by percentile
    data = np.squeeze(data[cmap])
    cmap = cmap[np.all([data > np.percentile(data,min_percentile), data < np.percentile(data,max_percentile)], axis=0)]     
    distance = distance[np.all([data > np.percentile(data,min_percentile), data < np.percentile(data,max_percentile)], axis=0)]     
    
    return cmap

# get cmap intersection
cmap = []
cmap_hemi = []    
for i in range(len(hemi)):
    cmap_multipol  = []
    for j in multipol:
        cmap_layer = []
            
        # load data
        input_cmap = os.path.join(path_ortho,hemi[i]+"."+file_patch+".patch.flat.cmap.nii")
        input_distance = os.path.join(path_ortho,hemi[i]+"."+file_patch+".iso_distance.nii")
        input_data = os.path.join(path_data,hemi[i]+".c_multipol_"+str(j)+"_def-img_layer"+str(layer_percentile)+"_def.mgh")
                
        c = get_percentile(input_data, input_cmap, input_distance, min_percentile, max_percentile, min_distance, max_distance)
        cmap_layer.append(c)
    
        temp = cmap_layer[0].copy()
        for n in range(len(cmap_layer)-1):
            temp = list(set(temp).intersection(cmap_layer[n+1]))
            
        cmap_multipol.append(temp)
    
    cmap_hemi = cmap_multipol[0].copy()
    for n in range(len(cmap_multipol)-1):
        cmap_hemi = list(set(cmap_hemi).intersection(cmap_multipol[n+1]))
    
    cmap.append(cmap_hemi)
    
# get distances
cmap_temp = []
distance_temp = []
distance = []
for i in range(len(hemi)):
    
    # load data
    input_cmap = os.path.join(path_ortho,hemi[i]+"."+file_patch+".patch.flat.cmap.nii")
    input_distance = os.path.join(path_ortho,hemi[i]+"."+file_patch+".iso_distance.nii")
    
    c = nb.load(input_cmap).get_fdata().astype(int)
    d = nb.load(input_distance).get_fdata()
    cmap_temp = np.append(cmap_temp, c)
    distance_temp = np.append(distance_temp, d)

    cmap_temp, ind = np.unique(cmap_temp, return_index=True)
    distance_temp = distance_temp[ind]

    distance_temp2 = np.zeros(len(cmap[i]))
    for j in range(len(cmap[i])):
        distance_temp2[j] = distance_temp[np.where(cmap_temp == cmap[i][j])[0][0]]
    
    distance.append(distance_temp2)

# get MTF data for each layer
for l in layer:

    # get multi data and corresponding frequency
    M = np.zeros((len(cmap[0])+len(cmap[1]),len(multipol)))
    F = np.zeros((len(cmap[0])+len(cmap[1]),len(multipol)))
    
    M = []
    F = []
    for i in multipol:
        for j in range(len(hemi)):
            input_data = os.path.join(path_data,hemi[j]+".c_multipol_"+str(i)+"_def-img_layer"+str(l)+"_def.mgh")
            data = nb.load(input_data).get_fdata()
            data = np.squeeze(data[cmap[j]])
            
            # get cortical frequency
            cortical_frequency = i / (2 * distance[j])
            
            # get final data M(F)
            M = np.append(M, data)
            F = np.append(F, cortical_frequency)

    # bin data
    bins = np.linspace(0,np.max(F),bin_size+1)
    bins = bins[0:-1]
    F_bin = np.digitize(F,bins) 
    
    M_bin_mean = np.full(len(bins), np.nan) 
    M_bin_std = np.full(len(bins), np.nan) 
    F_bin_mean = np.full(len(bins), np.nan) 
    F_bin_std = np.full(len(bins), np.nan) 
    for i in range(len(bins)):
        if np.size(M[F_bin == i+1]) > 0:
            M_bin_mean[i] = np.nanmean(M[F_bin == i+1])
            M_bin_std[i] = np.nanstd(M[F_bin == i+1])
            F_bin_mean[i] = np.nanmean(F[F_bin == i+1])
            F_bin_std[i] = np.nanstd(F[F_bin == i+1])
    
    # remove nans
    M_bin_mean = M_bin_mean[~np.isnan(M_bin_mean)]
    M_bin_std = M_bin_std[~np.isnan(M_bin_std)]
    F_bin_mean = F_bin_mean[~np.isnan(F_bin_mean)]
    F_bin_std = F_bin_std[~np.isnan(F_bin_std)]
    
    # gaussian fit
    mean = sum(F_bin_mean * M_bin_mean) / len(F_bin_mean)
    sigma = np.sqrt(sum(M_bin_mean * (F_bin_mean - mean)**2) / len(F_bin_mean))
    c = 0 
    
    # Define a gaussian function with offset
    def gauss(x, beta, sigma, c):
        return beta * np.exp(-x**2 / (2*sigma**2)) + c
    
    popt,pcov = curve_fit(gauss,F_bin_mean,M_bin_mean,p0=[1,sigma,c])
    
    # plot gaussian fit
    plt.errorbar(F_bin_mean,M_bin_mean,M_bin_std,F_bin_std,'b+:',label='data')
    plt.plot(F_bin_mean,gauss(F_bin_mean,*popt),'ro:',label='fit')
    plt.legend()
    plt.xlabel('Cortical frequency in cycles/mm')
    plt.ylabel('BOLD in a.u.')
    plt.show()

    # save data
    F_out = np.linspace(0,np.max(F_bin_mean),1000)
    multi_out = gauss(F_out,*popt)
    np.save(os.path.join(path_output,"x_"+str(l)),F_out)
    np.save(os.path.join(path_output,"y_"+str(l)),multi_out)
    np.save(os.path.join(path_output,"sigma_"+str(l)),popt[1])