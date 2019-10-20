"""
Merge outliers from vaso

Regressors of no interest from spatial correlation analysis and motion correction are derived from
single contrasts separtely. Here, outlier lists are merged together.

created by Daniel Haenelt
Date created: 27-09-2019             
Last modified: 06-10-2019  
"""
import sys
import os
import numpy as np

input_outlier_vaso = [
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_1/logfiles/outlier_regressor_vaso.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_2/logfiles/outlier_regressor_vaso.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_3/logfiles/outlier_regressor_vaso.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_4/logfiles/outlier_regressor_vaso.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_5/logfiles/outlier_regressor_vaso.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_6/logfiles/outlier_regressor_vaso.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_7/logfiles/outlier_regressor_vaso.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_8/logfiles/outlier_regressor_vaso.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_9/logfiles/outlier_regressor_vaso.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_10/logfiles/outlier_regressor_vaso.txt",
        ]
input_outlier_bold = [
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_1/logfiles/outlier_regressor_bold.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_2/logfiles/outlier_regressor_bold.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_3/logfiles/outlier_regressor_bold.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_4/logfiles/outlier_regressor_bold.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_5/logfiles/outlier_regressor_bold.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_6/logfiles/outlier_regressor_bold.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_7/logfiles/outlier_regressor_bold.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_8/logfiles/outlier_regressor_bold.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_9/logfiles/outlier_regressor_bold.txt",
        "/data/pt_01880/Experiment1_ODC/p3/odc/VASO3/Run_10/logfiles/outlier_regressor_bold.txt",
        ]

""" do not edit below """

for i in range(len(input_outlier_vaso)):
    
    # set output path and filename
    path_output = os.path.dirname(input_outlier_vaso[i])
    name_output = os.path.basename(input_outlier_vaso[i]).split("_")[0]+"_regressor_merge.txt"
    
    # load outlier
    outlier_vaso = np.loadtxt(input_outlier_vaso[i]).astype(int)
    outlier_bold = np.loadtxt(input_outlier_bold[i]).astype(int)
    
    # double entries
    outlier_vaso = np.repeat(outlier_vaso,2)
    outlier_bold = np.repeat(outlier_bold,2)

    # check same length of both outlier files
    if len(outlier_vaso) != len(outlier_bold):
        sys.exit("Outlier files do not have the same length in run "+str(i)+"!")

    # merge vaso and bold outliers
    outlier_merge = outlier_vaso + outlier_bold
    outlier_merge[outlier_merge > 1] = 1

    # save merged regressor
    np.savetxt(os.path.join(path_output,name_output),outlier_merge,'%i')