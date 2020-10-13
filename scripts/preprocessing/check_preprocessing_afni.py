# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys
import shutil as sh

# external inputs
import numpy as np
import matplotlib.pyplot as plt

# local inputs
from fmri_tools.io.get_filename import get_filename
from fmri_tools.utils.get_mean import get_mean
from fmri_tools.registration.clean_ana import clean_ana
from fmri_tools.registration.mask_ana import mask_ana
from fmri_tools.registration.mask_epi import mask_epi


"""
Outlier detection

This script looks for outliers in a functional time series using the Afni 
function 3dToutcount. The fraction of outliers found within a volume defined by 
a mask are written out for each time point. Additionally, a graphical 
visualization, a regressor of no interest and a short summary are written.

The script needs an installation of fsl and afni.

created by Daniel Haenelt
Date created: 19-02-2020 
Last modified: 13-10-2020  
"""

input_epi = [
    "/data/pt_01880/Experiment1_ODC/p4/retinotopy3/pol_anticlock/uadata.nii",
    ]
input_t1 = [
    "/data/pt_01880/Experiment1_ODC/p4/anatomy/S7_MP2RAGE_0p7_T1_Images_2.45.nii",
    ]
input_mask = [
    "/data/pt_01880/Experiment1_ODC/p4/anatomy/skull/skullstrip_mask_enhanced.nii",
    ]

# parameters
niter_mask = 3
sigma_mask = 3
qthr = 0.001
pthr = 0.01
cleanup = False

# do not edit below

# check length of input arrays
if len(input_epi) != len(input_t1) or len(input_t1) != len(input_mask):
    sys.exit("The lengths of input arrays are not consistent!")

for i in range(len(input_epi)):
    
    # get fileparts from first input entry
    path_epi, name_epi, ext_epi = get_filename(input_epi[i])
    _, _, ext_t1 = get_filename(input_t1[i])
    _, _, ext_mask = get_filename(input_mask[i])
    
    # make folder structure
    path_output = os.path.join(path_epi, "outlier_afni")
    path_temp = os.path.join(path_output, "temp")
    
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    if not os.path.exists(path_temp):
        os.makedirs(path_temp)

    # copy input files into temporary folder
    sh.copyfile(input_t1[i], os.path.join(path_temp,"T1"+ext_t1))
    sh.copyfile(input_mask[i], os.path.join(path_temp,"mask"+ext_mask))

    # get mean time series
    get_mean(input_epi[i], path_temp, "epi", type="mean")

    # new filenames
    file_epi_mean = os.path.join(path_temp, "mean_epi"+ext_epi)
    file_t1 = os.path.join(path_temp,"T1"+ext_t1)
    file_mask = os.path.join(path_temp,"mask"+ext_mask)

    # get mask
    clean_ana(file_t1, 1000.0, 4095.0, overwrite=True) # clean ana
    mask_ana(file_t1, file_mask, background_bright=False) # mask t1
    mask_epi(file_epi_mean, os.path.join(path_temp,"pT1"+ext_t1), file_mask, niter_mask, sigma_mask)
    
    # get outlier count within mask
    os.system("3dToutcount " + \
          "-mask " + os.path.join(path_temp,"mask_def-img3.nii.gz") + " " + \
          "-fraction " + \
          "-qthr " +str(qthr) + " " + input_epi[i] + " " + \
          " > " + os.path.join(path_output,"outlier_afni.txt"))
    
    # make plot
    log_data = np.loadtxt(os.path.join(path_output,"outlier_afni.txt"))

    # plots
    fig, ax = plt.subplots()
    ax.plot(log_data, "r")
    ax.set_xlabel("Time in TR")
    ax.set_ylabel("Outlier count")
    ax.set_title("Fraction of outliers found at each time step")
    ax.hlines(pthr,0,len(log_data),linestyle="dashed")
    fig.savefig(os.path.join(path_output,"outlier_afni.png"), format='png', bbox_inches='tight')

    # make regressor of no interest
    log_outlier = np.zeros_like(log_data)
    log_outlier[log_data > pthr] = 1
    np.savetxt(os.path.join(path_output,"outlier_afni_regressor_"+name_epi+".txt"), log_outlier, fmt="%d")

    # make summary textfile
    file = open(os.path.join(path_output,"outlier_afni_summary.txt"),"w")
    file.write("qthr: "+str(qthr)+"\n")
    file.write("pthr: "+str(pthr)+"\n")
    file.write("Number of outliers: "+str(np.sum(log_outlier))) 
    file.close()

    # clean intermediate files
    if cleanup:
        sh.rmtree(path_temp, ignore_errors=True)
