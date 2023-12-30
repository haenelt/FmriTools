# -*- coding: utf-8 -*-
"""
Gradient nonlinearity correction

This scripts calls the HCP toolbox to correct for gradient nonlinearities in the
input volume. Further explanations can be found in a separate readme. The script
needs an installation of fsl.

Dependencies:
- gradunwarp: https://github.com/Washington-University/gradunwarp
- please see the readme in the current folder for installation instructions

"""

import os

import fmri_tools

from ..io.filename import get_filename
from ..preprocessing.gnl_correction import gnl_correction

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
file_bash = os.path.join(
    os.path.dirname(fmri_tools.__file__), "preprocessing", "apply_grad.sh"
)

for i in range(len(file_in)):
    # get filename
    path_output, _, _ = get_filename(file_in[i])

    # gnl correction
    gnl_correction(
        file_in[i],
        file_bash,
        file_coeff,
        python3_env,
        python2_env,
        path_output,
        cleanup,
    )
