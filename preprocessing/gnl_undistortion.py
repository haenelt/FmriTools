"""
Gradient nonlinearity correction

This scripts calls the HCP toolbox to correct for gradient nonlinearities in the input volume.

Before running the script, login to queen via ssh and set the fsl environment by calling FSL in 
the terminal.

created by Daniel Haenelt
Date created: 10-01-2020
Last modified: 10-01-2020
"""
from lib.preprocessing.gnl_correction import gnl_correction

# input
input = [
    "/data/pt_01880/odc_temp/gnl/S23_MP2RAGE_0p7_UNI_Images_2.45.nii",
    ]

path_output = [
    "/data/pt_01880/odc_temp/gnl",
    ]

file_bash = "/home/raid2/haenelt/projects/gradunwarp/apply_grad.sh"
file_coeff = "/home/raid2/haenelt/projects/gradunwarp/7t_coeff.grad"
python3_env = "daniel"
python2_env = "daniel2"
cleanup = True

""" do not edit below """

for i in range(len(input)):
    gnl_correction(input[i], file_bash, file_coeff, python3_env, python2_env, path_output[i], 
                   cleanup)