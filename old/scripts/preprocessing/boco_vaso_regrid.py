# -*- coding: utf-8 -*-
"""
Bold correction of VASO data

This scripts corrects a vaso time series for bold contamination. First, both
time series are upsampled to a common time grid. BOLD correction is performed by
dividing both time series. In the end, unrealistic vaso values are removed. The
script needs an installation of afni.

"""

import os

import nibabel as nb
import numpy as np

from ..io.filename import get_filename
from ..utils.regrid_time_series import regrid_time_series

# input data
img_vaso = [
    "/data/pt_01880/temp/try2/Run_1/uvaso.nii",
    "/data/pt_01880/temp/try2/Run_2/uvaso.nii",
    "/data/pt_01880/temp/try2/Run_3/uvaso.nii",
    "/data/pt_01880/temp/try2/Run_4/uvaso.nii",
    "/data/pt_01880/temp/try2/Run_5/uvaso.nii",
    "/data/pt_01880/temp/try2/Run_6/uvaso.nii",
    "/data/pt_01880/temp/try2/Run_7/uvaso.nii",
    "/data/pt_01880/temp/try2/Run_8/uvaso.nii",
    "/data/pt_01880/temp/try2/Run_9/uvaso.nii",
    "/data/pt_01880/temp/try2/Run_10/uvaso.nii",
]

img_bold = [
    "/data/pt_01880/temp/try2/Run_1/ubold.nii",
    "/data/pt_01880/temp/try2/Run_2/ubold.nii",
    "/data/pt_01880/temp/try2/Run_3/ubold.nii",
    "/data/pt_01880/temp/try2/Run_4/ubold.nii",
    "/data/pt_01880/temp/try2/Run_5/ubold.nii",
    "/data/pt_01880/temp/try2/Run_6/ubold.nii",
    "/data/pt_01880/temp/try2/Run_7/ubold.nii",
    "/data/pt_01880/temp/try2/Run_8/ubold.nii",
    "/data/pt_01880/temp/try2/Run_9/ubold.nii",
    "/data/pt_01880/temp/try2/Run_10/ubold.nii",
]

# parameters
TR_old = 5  # effective TR of bold+vaso
TR_new = 3  # TR of upsampled corrected time series
start_bold = 0  # start of bold block in s
start_vaso = 0  # start of vaso block in s
vaso_threshold = 6  # threshold unrealistic intensities
nvol_remove = 0  # number of volumes removed at the end of the time series

# do not edit below

for bold, vaso in zip(img_bold, img_vaso):
    # get filenames
    path_bold, name_bold, ext_bold = get_filename(bold)
    path_vaso, name_vaso, ext_vaso = get_filename(vaso)

    # upsample time series
    regrid_time_series(
        bold, path_bold, TR_old, TR_new, t_start=start_bold, nvol_remove=nvol_remove
    )
    regrid_time_series(
        vaso, path_vaso, TR_old, TR_new, t_start=start_vaso, nvol_remove=nvol_remove
    )

    # new filenames
    file_bold = os.path.join(path_bold, name_bold + "_upsampled" + ext_bold)
    file_vaso = os.path.join(path_vaso, name_vaso + "_upsampled" + ext_vaso)
    file_vaso_corrected = os.path.join(
        path_vaso, name_vaso + "_upsampled_corrected" + ext_vaso
    )

    # load bold data
    bold_data = nb.load(file_bold)
    bold_array = bold_data.get_fdata()

    # load vaso data
    vaso_data = nb.load(file_vaso)
    vaso_array = vaso_data.get_fdata()

    # bold correction
    vaso_array = np.divide(vaso_array, bold_array)

    # remove nans and infs
    vaso_array[np.isnan(vaso_array)] = 0
    vaso_array[np.isinf(vaso_array)] = 0

    # clean vaso data that are unrealistic
    vaso_array[vaso_array < 0] = 0
    vaso_array[vaso_array >= vaso_threshold] = vaso_threshold

    # write output
    output = nb.Nifti1Image(vaso_array, vaso_data.affine, vaso_data.header)
    nb.save(output, file_vaso_corrected)
