"""
Merge outliers from vaso

Regressors of no interest from spatial correlation analysis and motion correction are derived from
single contrasts separtely. Here, outlier lists are merged together.

created by Daniel Haenelt
Date created: 27-09-2019             
Last modified: 27-09-2019  
"""
import sys
import os
import numpy as np

input_outlier_vaso = ["/data/pt_01880/Experiment1_ODC/p1/odc/SE_EPI2/Run_1/logfiles/correlation_regressor_vaso.txt"]
input_outlier_bold = ["/data/pt_01880/Experiment1_ODC/p1/odc/SE_EPI2/Run_1/logfiles/correlation_regressor_bold.txt"]


""" do not edit below """

for i in range(len(input_outlier_vaso)):
    
    # set output path and filename
    path_output = os.path.dirname(input_outlier_vaso[i])
    name_output = os.path.basename(input_outlier_vaso[i]).split("_")[0]+"_regressor_merge.txt"
    
    # load outlier
    outlier_vaso = np.loadtxt(input_outlier_vaso[i]).astype(int)
    outlier_bold = np.loadtxt(input_outlier_bold[i]).astype(int)

    # check same length of both outlier files
    if len(outlier_vaso) != len(outlier_bold):
        sys.exit("Outlier files do not have the same length in run "+str(i)+"!")

    outlier_merge = []
    for j in range(len(outlier_vaso)):
        outlier_merge = np.append(outlier_merge, [outlier_vaso[j], outlier_bold[j]]).astype(int)

    # save merged regressor
    np.savetxt(os.path.join(path_output,name_output),outlier_merge,'%i')
