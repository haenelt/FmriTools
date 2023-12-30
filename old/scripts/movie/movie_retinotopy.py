# -*- coding: utf-8 -*-
"""
Make retinotopy movie

The purpose of the following script is to make a movie to illustrate the
travelling waves on a flattened patch.

"""

import glob
import os

import numpy as np

from ..img.get_gif import get_gif
from ..img.get_retinotopy_images import get_retinotopy_images

# input files
input_patch = "/data/pt_01880/V2STRIPES/p6/anatomy/dense/lh.occip3.patch.flat"
input_vfs = "/data/pt_01880/lh.fieldsign.mgh"
input_phase = (
    "/data/pt_01880/V2STRIPES/p6/retinotopy/avg/surf/lh.pol_phase_avg_def_layer5.mgh"
)
input_snr = (
    "/data/pt_01880/V2STRIPES/p6/retinotopy/avg/surf/lh.pol_snr_avg_def_layer5.mgh"
)
input_white = "/data/pt_01880/V2STRIPES/p6/anatomy/dense/lh.white"

# retinotopy img input
hemi = "lh"
path_output = "/data/pt_01880"

# gif parameters
name_output = "pol"
nsteps = 1
duration = 0.025

# do not edit below


def _natsort(file_list):
    """Sort a list of file names based on their first natural number in the file name.
    This sorting is equivalent to the natsorted() function in the natsort package.

    Parameters
    ----------
    file_list : list
        List of file names.

    Returns
    -------
    list
        Sorted list of file names.
    """
    file_number = []
    for file_ in file_list:
        file_number.append([int(i) for i in file_.split() if i.isdigit()][0])
    file_ind = np.argsort(file_number)
    return [file_list[i] for i in file_ind]


# get single frames of the travelling wave
get_retinotopy_images(
    input_patch,
    input_vfs,
    input_phase,
    input_snr,
    input_white,
    hemi,
    path_output,
    img_res=0.2,
    theta=0,
    alpha=2,
    buffer=0,
    phase_fwhm=4,
    sigma=50,
    cleanup=False,
)

# make gif
img_files = glob.glob(os.path.join(path_output, "img", "*"))
img_files = _natsort(img_files)
get_gif(img_files, path_output, name_output, nsteps, duration)
