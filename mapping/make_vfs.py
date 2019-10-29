"""
Make visual fieldsign map

The purpose of the following script is to make a visual fieldsign map from sampled phase-encoded
retinotopy data.

Before running the script, login to queen via ssh and set the freesurfer environments by calling 
FREESURFER in the terminal.

created by Daniel Haenelt
Date created: 14-02-2019
Last modified: 29-08-2019
"""
import os
from lib.mapping.get_vfs import get_vfs

# input files
hemi = "lh"
input_sphere = "/data/pt_01880/Experiment1_ODC/p3/anatomy/dense/lh.sphere"
input_white = "/data/pt_01880/Experiment1_ODC/p3/anatomy/dense/lh.white"
input_patch = "/data/pt_01880/Experiment1_ODC/p3/anatomy/dense/lh.occip.patch.flat"
input_aparc = "/data/pt_01880/Experiment1_ODC/p3/anatomy/freesurfer/label/lh.aparc.annot"
ecc_real = "/data/pt_01880/Experiment1_ODC/p3/retinotopy/avg/surf/lh.ecc_real_avg_def-img_layer5_def.mgh"
ecc_imag = "/data/pt_01880/Experiment1_ODC/p3/retinotopy/avg/surf/lh.ecc_imag_avg_def-img_layer5_def.mgh"
pol_real = "/data/pt_01880/Experiment1_ODC/p3/retinotopy/avg/surf/lh.pol_real_avg_def-img_layer5_def.mgh"
pol_imag = "/data/pt_01880/Experiment1_ODC/p3/retinotopy/avg/surf/lh.pol_imag_avg_def-img_layer5_def.mgh"
path_output = "/data/pt_01880/Experiment1_ODC/p3/retinotopy/avg/surf/"

""" do not edit below """

get_vfs(input_sphere, input_white, input_patch, input_aparc, hemi, ecc_real, ecc_imag, pol_real, 
            pol_imag, path_output, fwhm_ecc = 4.0, fwhm_pol = 2.0, fwhm_vfs = 8.0, cleanup=True)

# remove generated tmp.dat file
os.remove(os.path.join(os.getcwd(),"tmp.dat"))
