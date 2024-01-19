# -*- coding: utf-8 -*-
"""
Temporal signal-to-noise ratio of functional time series

This scripts calculates temporal signal-to-noise ratio of a set of time series.
TSNR os single time series and the overall mean are written.

"""

import datetime
import os

import nibabel as nb
import numpy as np

from ..io.filename import get_filename
from ..matlab import MatlabCommand
from ..utils.metrics import calc_tsnr

# input data
img_input = [
    "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI1/Run_1/udata.nii",
    "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI1/Run_2/udata.nii",
    "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI1/Run_3/udata.nii",
    "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI1/Run_4/udata.nii",
    "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI1/Run_5/udata.nii",
    "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI1/Run_6/udata.nii",
    "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI1/Run_7/udata.nii",
    "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI1/Run_8/udata.nii",
    "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI1/Run_9/udata.nii",
    "/data/pt_01880/Experiment1_ODC/p4/odc/GE_EPI1/Run_10/udata.nii",
]

# parameters
tsnr_threshold = 200  # threshold unrealistic high tsnr values
use_highpass = False
TR = 3  # repetition time in s
cutoff_highpass = 270  # cutoff in s for baseline correction
name_sess = "GE_EPI1"
name_output = ""

# do not edit below

# get path from first entry
path_file, _, _ = get_filename(img_input[0])

# make output folder
path_output = os.path.join(
    os.path.dirname(os.path.dirname(path_file)), "results", "tsnr", "native"
)
if not os.path.exists(path_output):
    os.makedirs(path_output)

# get image header information
data_img = nb.load(img_input[0])
data_img.header["dim"][0] = 3
data_img.header["dim"][4] = 1
header = data_img.header
affine = data_img.affine

# get image dimension
dim = data_img.header["dim"][1:4]

mean_tsnr = np.zeros(dim)
for i in range(len(img_input)):
    # get filename
    path_file, name_file, ext_file = get_filename(img_input[i])

    # highpass filter time series
    if use_highpass:
        matlab = MatlabCommand(
            "ft_baseline_correction",
            os.path.join(path_file, name_file + ext_file),
            TR,
            cutoff_highpass,
        )
        matlab.run()

        # change input to highpass filtered time series
        name_file = "b" + name_file

    # sum volumes for each run
    mean_tsnr += calc_tsnr(
        os.path.join(path_file, name_file + ext_file),
        tsnr_max=tsnr_threshold,
        write_output=True,
        path_output=path_file,
    )

# divide by number of runs
mean_tsnr /= len(img_input)

# name of output files
if len(name_output) and len(name_sess):
    fileOUT = os.path.join(
        path_output, "tsnr_" + name_output + "_" + name_sess + ".nii"
    )
elif len(name_output) and not len(name_sess):
    fileOUT = os.path.join(path_output, "tsnr_" + name_output + ".nii")
elif not len(name_output) and len(name_sess):
    fileOUT = os.path.join(path_output, "tsnr_" + name_sess + ".nii")
else:
    fileOUT = os.path.join(path_output, "tsnr.nii")

# write output
output = nb.Nifti1Image(mean_tsnr, affine, header)
nb.save(output, fileOUT)

# write log
fileID = open(os.path.join(path_output, "tsnr_info.txt"), "a")
fileID.write(
    "script executed: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n"
)
fileID.write("session: " + name_sess + "\n")
fileID.write("basename: " + name_output + "\n")
fileID.write("tsnr threshold: " + str(tsnr_threshold) + "\n")
fileID.write("TR: " + str(TR) + "\n")
fileID.write("highpass: " + str(use_highpass) + "\n")
fileID.write("cutoff highpass: " + str(cutoff_highpass) + "\n")
fileID.close()
