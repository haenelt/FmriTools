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
input = ["/data/pt_01880/Experiment5_Mono/p1/flicker/GE_EPI1/Run_1/data.nii",
         "/data/pt_01880/Experiment5_Mono/p1/flicker/GE_EPI1/Run_2/data.nii",
         "/data/pt_01880/Experiment5_Mono/p1/flicker/GE_EPI1/Run_3/data.nii",
         "/data/pt_01880/Experiment5_Mono/p1/flicker/GE_EPI1/Run_4/data.nii",
         "/data/pt_01880/Experiment5_Mono/p1/flicker/GE_EPI1/Run_5/data.nii",
         "/data/pt_01880/Experiment5_Mono/p1/flicker/GE_EPI2/Run_1/data.nii",
         "/data/pt_01880/Experiment5_Mono/p1/flicker/GE_EPI2/Run_2/data.nii",
         "/data/pt_01880/Experiment5_Mono/p1/flicker/GE_EPI2/Run_3/data.nii",
         "/data/pt_01880/Experiment5_Mono/p1/flicker/GE_EPI2/Run_4/data.nii",
         "/data/pt_01880/Experiment5_Mono/p1/flicker/GE_EPI2/Run_5/data.nii",
         "/data/pt_01880/Experiment5_Mono/p1/flicker/GE_EPI2/Run_6/data.nii",
         "/data/pt_01880/Experiment5_Mono/p1/flicker/GE_EPI2/Run_7/data.nii",
         ]
TR_old = 3
TR_new = 1
order = "interleaved"

""" do not edit below """

for i in range(len(input)):
    slice_timing_correction(input[i], TR_old, TR_new, order, prefix="a")
