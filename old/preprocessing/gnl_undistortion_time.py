# -*- coding: utf-8 -*-
"""
Gradient nonlinearity correction of time series

This scripts calls the HCP toolbox to correct for gradient nonlinearities in the
input time series. The corrected time series gets the suffix _gnlcorr. The
script needs an installation of fsl.

Dependencies:
- gradunwarp: https://github.com/Washington-University/gradunwarp
- please see the readme in the current folder for installation instructions

"""

import os
import shutil as sh

import nibabel as nb
import numpy as np

import fmri_tools

from ..io.filename import get_filename
from ..preprocessing.gradient import gnl_correction
from ..registration.fsl import apply_warp
from ..utils.roi import extract_vol

# input
file_in = [
    "/data/pt_01880/Experiment3_Stripes/p3/colour/GE_EPI1/Run_1/data.nii",
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

for file_ in file_in:
    # get fileparts of input
    path_file, name_file, ext_file = get_filename(file_)

    # filenames
    file_vol0 = os.path.join(path_file, name_file + "_vol0" + ext_file)
    file_out = os.path.join(path_file, name_file + "_gnlcorr" + ext_file)

    # extract first volume
    extract_vol(file_, file_vol0, 0, 1)

    # exexute gnl correction
    gnl_correction(
        file_vol0, file_bash, file_coeff, python3_env, python2_env, path_file, False
    )

    # apply warp to first volume
    apply_warp(file_, os.path.join(path_file, "grad", "warp.nii.gz"), file_out)

    # normalise warped output image to initial intensity range
    data_img = nb.load(file_)
    data_array = data_img.get_fdata()
    max_data = np.max(data_array)
    min_data = np.min(data_array)

    data_img = nb.load(file_out)
    data_array = data_img.get_fdata()
    data_array[data_array < min_data] = 0
    data_array[data_array > max_data] = max_data

    output = nb.Nifti1Image(data_array, data_img.affine, data_img.header)
    nb.save(output, file_out)

    # remove vol0
    os.remove(file_vol0)
    os.remove(os.path.join(path_file, name_file + "_vol0_gnlcorr" + ext_file))

    # clean intermediate files
    if cleanup:
        sh.rmtree(os.path.join(path_file, "grad"), ignore_errors=True)
