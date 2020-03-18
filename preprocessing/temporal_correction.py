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
input = ["/data/pt_01880/temp_slice/slice_yes/Run_1/data.nii",
         "/data/pt_01880/temp_slice/slice_yes/Run_2/data.nii",
         "/data/pt_01880/temp_slice/slice_yes/Run_3/data.nii",
         "/data/pt_01880/temp_slice/slice_yes/Run_4/data.nii",
         "/data/pt_01880/temp_slice/slice_yes/Run_5/data.nii",
         "/data/pt_01880/temp_slice/slice_yes/Run_6/data.nii",
         "/data/pt_01880/temp_slice/slice_yes/Run_7/data.nii",
         "/data/pt_01880/temp_slice/slice_yes/Run_8/data.nii",
         "/data/pt_01880/temp_slice/slice_yes/Run_9/data.nii",
         "/data/pt_01880/temp_slice/slice_yes/Run_10/data.nii",
         ]
TR_old = 3
TR_new = 2
order = "interleaved"

""" do not edit below """

for i in range(len(input)):
    slice_timing_correction(input[i], TR_old, TR_new, order, prefix="a")
