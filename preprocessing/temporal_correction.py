"""
Temporal correction

This scripts performs slice timing correction and temporal regridding on a nifti time series.

Before running the script, login to queen via ssh and set the afni environment by calling AFNI in 
the terminal.

created by Daniel Haenelt
Date created: 11-03-2020
Last modified: 11-03-2020
"""
from lib.preprocessing.slice_timing_correction import slice_timing_correction

# input
input = ["/data/pt_01880/Experiment5_Mono/p1/flicker/GE_EPI1/Run_1/data.nii"]
TR_old = [3]
order = ["interleaved"]
TR_new = 1

""" do not edit below """

for i in range(len(input)):
    slice_timing_correction(input[i], TR_old[i], TR_new, order[i], prefix="a")