"""
Make visual fieldsign map

The prupose of the following script is to make a movie to illustrate the travelling waves on a 
flattened patch.

created by Daniel Haenelt
Date created: 14-02-2019
Last modified: 18-02-2019
"""
from lib.mapping.get_vfs import get_vfs

# input files
input_sphere = "/data/pt_01880/V2STRIPES/p6/anatomy/dense/rh.sphere"
input_white = "/data/pt_01880/V2STRIPES/p6/anatomy/dense/rh.white"
input_patch = "/data/pt_01880/V2STRIPES/p6/anatomy/dense/rh.occip3.patch.flat"
input_aparc = "/data/pt_01880/V2STRIPES/p6/anatomy/freesurfer/label/rh.aparc.annot"
hemi = "rh"
ecc_real = "/data/pt_01880/V2STRIPES/p6/retinotopy/avg/surf/rh.ecc_real_avg_def_layer5.mgh"
ecc_imag = "/data/pt_01880/V2STRIPES/p6/retinotopy/avg/surf/rh.ecc_imag_avg_def_layer5.mgh"
pol_real = "/data/pt_01880/V2STRIPES/p6/retinotopy/avg/surf/rh.pol_real_avg_def_layer5.mgh"
pol_imag = "/data/pt_01880/V2STRIPES/p6/retinotopy/avg/surf/rh.pol_imag_avg_def_layer5.mgh"
path_output = "/data/pt_01880"

""" do not edit below """

get_vfs(input_sphere, input_white, input_patch, input_aparc, hemi, ecc_real, ecc_imag, pol_real, 
            pol_imag, path_output, fwhm_ecc = 4.0, fwhm_pol = 2.0, fwhm_vfs = 8.0, cleanup=True)