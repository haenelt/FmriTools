# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np

# local inputs
from fmri_tools.io.get_filename import get_filename


"""
Upsample vaso outliers

Since the corrected vaso time series is regridded to a new TR, outlier lists 
e.g. to censor single volumes based on the realignment procedure, are 
transformed to the corresponding timepoints in the resampled timeseries.

created by Daniel Haenelt
Date created: 26-05-2020           
Last modified: 12-10-2020  
"""

# input
file_in = ["/data/pt_01880/Experiment1_ODC/p4/odc/VASO1/Run_3/outlier/outlier_regressor_merge.txt"]

# parameters
TR_old = 5 # effective TR of bold+vaso
TR_new = 3 # TR of upsampled bold corrected time series

# do not edit below

for i in range(len(file_in)):
    
    # set output path and filename
    path_output, _, _ = get_filename(file_in[i])
    name_output = "outlier_regressor_upsampled.txt"

    # load outlier textfile
    outlier = np.loadtxt(file_in[i])

    # merge bold and vaso outliers to one timepoint    
    outlier1 = outlier[::2]
    outlier2 = outlier[1::2]
    outlier_merge = outlier1 + outlier2
    outlier_merge[outlier_merge != 0] = 1

    # get time axes
    run_length = len(outlier_merge) * TR_old
    
    nt_old = int(run_length / TR_old)
    t_old = np.arange(0,nt_old) * TR_old
    
    nt_new = int(run_length / TR_new)
    t_new = np.arange(0,nt_new) * TR_new
    
    # find outliers
    n_outlier = np.where(outlier_merge == 1)[0]
    t_outlier = t_old[n_outlier]
    
    # get nearest time points in new time array
    outlier_new = np.zeros_like(t_new)
    for i in range(len(t_new)):
        for j in range(len(t_outlier)):
            if t_new[i] > t_outlier[j] - TR_new and t_new[i] < t_outlier[j] + TR_old:
                outlier_new[i] = 1

    # save merged regressor
    np.savetxt(os.path.join(path_output,name_output),outlier_new,'%i')
