# -*- coding: utf-8 -*-
"""
ALFF

The purpose of the following script is to compute ALFF and fALFF for a
resting-state time series. The time series is baseline corrected and nuisance
regressors are taken either from mean WM and CSF masks or based on biopac
recordings. Optionally, no nuisance regression is performed. The script needs an
installation of freesurfer, ants adn afni.

Dependencies (should be installed as SPM toolbox):
- PhysIO toolbox: https://www.nitrc.org/projects/physio/

"""

import os

import nibabel as nb

from ..matlab import MatlabCommand
from ..preprocessing.noise import nuisance_regressor
from ..processing.rest import get_alff
from ..segmentation.mask import mask_nuisance

# input
anatomy = ""  # T1w full brain anatomy (e.g. orig)
function = "/data/pt_01983/func/resting_state2/uadata.nii"  # baseline uncorrected
deformation = ""  # deformation ana -> epi
biopac_input = ""  # *.mat file
path_output = "/data/pt_01983/func/resting_state2/alff/native"

# parameters
TR = 3  # repetition time in s
cutoff_highpass = 270  # cutoff frequency for baseline correction in 1/Hz
nerode_wm = 1  # number of wm mask eroding iterations
nerode_csf = 1  # number of csf mask eroding iterations
hp_freq = 0.01  # highpass cutoff frequency (bandpass filter) in Hz
lp_freq = 0.08  # lowpass cutoff frequency (bandpass filter) in Hz

# analysis type
nuisance_regression = False
segmentation = False
biopac = False
cleanup = True

# do not edit below

# make output folder
if not os.path.exists(path_output):
    os.makedirs(path_output)

# get path and filenames
path = os.path.dirname(function)
file = os.path.basename(function)
bfile = "b" + file  # filename of baseline corrected time series
rfile = "r" + file  # filename of residual time series

# physiological noise regression
if nuisance_regression:
    # baseline correction
    matlab = MatlabCommand("ft_baseline_correction", function, TR, cutoff_highpass)
    matlab.run()

    if biopac:
        # get biopac regressors
        matlab = MatlabCommand(
            "ft_biopac_regressor",
            os.path.join(path, bfile),
            biopac_input,
            path_output,
            TR,
        )
        matlab.run()

    else:
        # get wm and csf mask
        mask_nuisance(
            anatomy,
            deformation,
            path_output,
            nerode_wm,
            nerode_csf,
            segmentation,
            cleanup,
        )

        # set mask to zero where function is equal to zero
        func_array = nb.load(function).get_fdata()
        func_array = func_array[:, :, :, 0]

        wm = nb.load(os.path.join(path_output, "wm_mask.nii.gz"))
        wm_array = wm.get_fdata()
        wm_array[func_array == 0] = 0
        output = nb.Nifti1Image(wm_array, wm.affine, wm.header)
        nb.save(output, os.path.join(path_output, "wm_mask.nii.gz"))

        csf = nb.load(os.path.join(path_output, "csf_mask.nii.gz"))
        csf_array = csf.get_fdata()
        csf_array[func_array == 0] = 0
        output = nb.Nifti1Image(csf_array, csf.affine, csf.header)
        nb.save(output, os.path.join(path_output, "csf_mask.nii.gz"))

        # get nuisance regressor
        nuisance_regressor(
            os.path.join(path, bfile),
            os.path.join(path_output, "wm_mask.nii.gz"),
            os.path.join(path_output, "csf_mask.nii.gz"),
            path_output,
        )

    # nuisance regression
    if cleanup:
        clean_glm = 1
    else:
        clean_glm = 0

    matlab = MatlabCommand(
        "ft_regress_physio",
        function,
        os.path.join(path_output, "nuisance_regressor.txt"),
        TR,
        cutoff_highpass,
        path_output,
        clean_glm,
    )
    matlab.run()

    # get alff
    get_alff(
        os.path.join(path_output, rfile), TR, path_output, hp_freq, lp_freq, cleanup
    )

    # remove baseline corrected time series
    if cleanup:
        os.remove(os.path.join(path, bfile))

else:
    # get alff
    get_alff(function, TR, path_output, hp_freq, lp_freq, cleanup)
