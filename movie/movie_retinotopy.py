"""
Make retinotopy movie

The purpose of the following script is to make a movie to illustrate the travelling waves on a 
flattened patch.

created by Daniel Haenelt
Date created: 14-02-2019
Last modified: 11-11-2019
"""
import os
import glob
import natsort
from lib.img.get_retinotopy_images import get_retinotopy_images
from lib.img.get_gif import get_gif

# input files
input_patch = "/data/pt_01880/V2STRIPES/p6/anatomy/dense/lh.occip3.patch.flat"
input_vfs = "/data/pt_01880/lh.fieldsign.mgh"
input_phase = "/data/pt_01880/V2STRIPES/p6/retinotopy/avg/surf/lh.pol_phase_avg_def_layer5.mgh"
input_snr = "/data/pt_01880/V2STRIPES/p6/retinotopy/avg/surf/lh.pol_snr_avg_def_layer5.mgh"
input_white = "/data/pt_01880/V2STRIPES/p6/anatomy/dense/lh.white"

# retinotopy img input
hemi = "lh"
path_output = "/data/pt_01880"

# gif parameters
name_output = "pol"
nsteps = 1
duration = 0.025

""" do not edit below """

# get single frames of the travelling wave
get_retinotopy_images(input_patch, input_vfs, input_phase, input_snr, input_white, hemi, 
                      path_output, img_res=0.2, theta=0, alpha=2, buffer=0, phase_fwhm=4, sigma=50,
                      cleanup=False)

# make gif
img_files = glob.glob(os.path.join(path_output,"img","*"))
img_files = natsort.natsorted(img_files)
get_gif(img_files, path_output, name_output, nsteps, duration)
