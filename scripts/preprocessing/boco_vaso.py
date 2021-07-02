# -*- coding: utf-8 -*-
"""
Bold correction of VASO data

This scripts corrects a vaso time series for bold contamination. First, both
time series are upsampled and the vaso time series is shifted by one time step.
BOLD correction is performed by dividing both time series. In the end,
unrealistic vaso values are removed. The script needs an installation of afni.

"""

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb

# input data
img_vaso = ["/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_1/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_2/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_3/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_4/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_5/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_6/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_7/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_8/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_9/uvaso.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_10/uvaso.nii",
            ]

img_bold = ["/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_1/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_2/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_3/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_4/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_5/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_6/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_7/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_8/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_9/ubold.nii",
            "/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_10/ubold.nii",
            ]

# parameters
TR = 2.5
vaso_threshold = 6

# do not edit below

for i in range(len(img_vaso)):
    # prepare path and filename
    path_vaso = os.path.dirname(img_vaso[i])
    file_vaso = os.path.splitext(os.path.basename(img_vaso[i]))[0]
    path_bold = os.path.dirname(img_bold[i])
    file_bold = os.path.splitext(os.path.basename(img_bold[i]))[0]

    # upsample vaso and bold time series
    os.system("3dUpsample -overwrite -datum short " +
              "-prefix " + os.path.join(path_vaso,
                                        file_vaso + "_upsampled.nii") +
              " -n 2 -input " + img_vaso[i])

    os.system("3dUpsample -overwrite -datum short " +
              "-prefix " + os.path.join(path_bold,
                                        file_bold + "_upsampled.nii") +
              " -n 2 -input " + img_bold[i])

    # load vaso data and shift in time
    vaso = nb.load(os.path.join(path_vaso, file_vaso + "_upsampled.nii"))
    vaso_array = vaso.get_fdata()
    vaso_array = vaso_array[:, :, :, :-1]
    vaso_array = np.concatenate((np.expand_dims(vaso_array[:, :, :, 0], axis=3),
                                 vaso_array), axis=3)

    # load bold data
    bold = nb.load(os.path.join(path_bold, file_bold + "_upsampled.nii"))
    bold_array = bold.get_fdata()

    # bold correction
    vaso_array = np.divide(vaso_array, bold_array)

    # remove nans and infs
    vaso_array[np.isnan(vaso_array)] = 0
    vaso_array[np.isinf(vaso_array)] = 0

    # clean vaso data that are unrealistic
    vaso_array[vaso_array < 0] = 0
    vaso_array[vaso_array >= vaso_threshold] = vaso_threshold

    output = nb.Nifti1Image(vaso_array, vaso.affine, vaso.header)
    nb.save(output, os.path.join(path_vaso,
                                 file_vaso + "_upsampled_corrected.nii"))

    # change TR in header
    os.system("3drefit " +
              "-TR " + str(TR) + " " +
              os.path.join(path_bold, file_bold + "_upsampled.nii"))

    os.system("3drefit " +
              "-TR " + str(TR) + " " +
              os.path.join(path_vaso, file_vaso + "_upsampled.nii"))

    os.system("3drefit " +
              "-TR " + str(TR) + " " +
              os.path.join(path_vaso, file_vaso + "_upsampled_corrected.nii"))
