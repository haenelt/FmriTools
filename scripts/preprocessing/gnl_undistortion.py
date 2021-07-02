# -*- coding: utf-8 -*-
"""
Gradient nonlinearity correction

This scripts calls the HCP toolbox to correct for gradient nonlinearities in the
input volume. Further explanations can be found in a separate readme. The script
needs an installation of fsl.

"""

# python standard library inputs
import os

# local inputs
import fmri_tools
from fmri_tools.preprocessing.gnl_correction import gnl_correction
from fmri_tools.io.get_filename import get_filename

# input
file_in = [
    "/data/pt_01880/Experiment3_Stripes/p2/anatomy/S5_MP2RAGE_0p7_INV1_2.45.nii",
    ]

file_coeff = "/data/hu_haenelt/projects/gradunwarp/7t_coeff.grad"
python3_env = "daniel"
python2_env = "daniel2"
cleanup = True

# do not edit below

# get path of bash file
file_bash = os.path.join(fmri_tools.__path__, 
                         "preprocessing", 
                         "apply_grad.sh")

for i in range(len(file_in)):
    
    # get filename
    path_output, _, _ = get_filename(file_in[i])
    
    # gnl correction
    gnl_correction(file_in[i],
                   file_bash,
                   file_coeff,
                   python3_env,
                   python2_env,
                   path_output,
                   cleanup)
