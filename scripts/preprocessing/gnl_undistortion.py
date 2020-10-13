# -*- coding: utf-8 -*-

# local inputs
from fmri_tools.preprocessing.gnl_correction import gnl_correction
from fmri_tools.io.get_filename import get_filename


"""
Gradient nonlinearity correction

This scripts calls the HCP toolbox to correct for gradient nonlinearities in the 
input volume.

The script needs an installation of fsl.

created by Daniel Haenelt
Date created: 10-01-2020
Last modified: 13-10-2020
"""

# input
input = [
    "/data/pt_01880/Experiment3_Stripes/p2/anatomy/S5_MP2RAGE_0p7_INV1_2.45.nii",
    ]

file_bash = "/data/hu_haenelt/projects/gradunwarp/apply_grad.sh"
file_coeff = "/data/hu_haenelt/projects/gradunwarp/7t_coeff.grad"
python3_env = "daniel"
python2_env = "daniel2"
cleanup = True

# do not edit below

for i in range(len(input)):
    
    # get filename
    path_output, _, _ = get_filename(input[i])
    
    # gnl correction
    gnl_correction(input[i], file_bash, file_coeff, python3_env, python2_env, path_output, cleanup)
