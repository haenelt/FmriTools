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
input_sphere = "/nobackup/eminem2/attar/3T_Connectom/pilot/subj-01/smri/mp2rage/2.0-analysis/5.0-surface/freesurfer/surf/lh.sphere"
input_white = "/nobackup/eminem2/attar/3T_Connectom/pilot/subj-01/smri/mp2rage/2.0-analysis/5.0-surface/freesurfer/surf/lh.white"
input_patch = "/nobackup/eminem2/attar/3T_Connectom/pilot/subj-01/smri/mp2rage/2.0-analysis/5.0-surface/freesurfer/surf/lh.occip.patch.flat"
input_aparc = "/nobackup/eminem2/attar/3T_Connectom/pilot/subj-01/smri/mp2rage/2.0-analysis/5.0-surface/freesurfer/label/lh.aparc.annot"
hemi = "lh"
ecc_real = "/nobackup/eminem2/attar/3T_Connectom/pilot/subj-01/fmri/retinotopy/avg/surf/lh.ecc_real_avg_def-img_mid_def.mgh"
ecc_imag = "/nobackup/eminem2/attar/3T_Connectom/pilot/subj-01/fmri/retinotopy/avg/surf/lh.ecc_imag_avg_def-img_mid_def.mgh"
pol_real = "/nobackup/eminem2/attar/3T_Connectom/pilot/subj-01/fmri/retinotopy/avg/surf/lh.pol_real_avg_def-img_mid_def.mgh"
pol_imag = "/nobackup/eminem2/attar/3T_Connectom/pilot/subj-01/fmri/retinotopy/avg/surf/lh.pol_imag_avg_def-img_mid_def.mgh"
path_output = "/nobackup/eminem2/attar/3T_Connectom/pilot/subj-01/fmri/retinotopy/avg/surf"

""" do not edit below """

get_vfs(input_sphere, input_white, input_patch, input_aparc, hemi, ecc_real, ecc_imag, pol_real, 
            pol_imag, path_output, fwhm_ecc = 4.0, fwhm_pol = 2.0, fwhm_vfs = 8.0, cleanup=True)