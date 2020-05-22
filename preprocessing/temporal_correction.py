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
input = ["/data/pt_01880/Experiment2_Rivalry/p1/rivalry/GE_EPI1/Run_1/data.nii",
         ]
TR_old = 1.5
TR_new = 1
order = "descending"

""" do not edit below """

for i in range(len(input)):
    slice_timing_correction(input[i], TR_old, TR_new, order, prefix="a")
