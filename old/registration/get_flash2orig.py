# -*- coding: utf-8 -*-

import os
import shutil as sh

from nighres.registration import apply_coordinate_mappings
from nipype.interfaces.fsl.preprocess import ApplyXFM

from ..io.vol import mri_convert
from ..registration.cmap import generate_coordinate_mapping
from ..registration.fsl import flirt
from ..registration.transform import scanner_transform


def get_flash2orig(file_flash, file_inv2, file_orig, path_output, cleanup=False):
    """Get FLASH to orig.

    This function computes the deformation field for the registration between a
    partial coverage GRE image and the freesurfer orig file. The following steps
    are performed: (1) set output folder structure, (2) get scanner transform
    GRE <-> inv2 and inv2 -> orig, (3) generate flash cmap, (4) apply scanner
    transform inv2 -> GRE, (5) get flirt registration GRE -> inv2, (6) apply
    flirt to GRE cmap, (7) apply scanner transform to GRE cmap, (8) apply final
    deformation to GRE. The function needs the FSL environment set.

    Parameters
    ----------
    file_flash : str
        Input path for GRE image.
    file_inv2 : str
        Input path for MP2RAGE INV2 image.
    file_orig : str
        Input path for freesurfer orig image.
    path_output : str
        Path where output is saved.
    cleanup : bool, optional
        Delete intermediate files. The default is False.

    Returns
    -------
    None.

    """

    # set folder structure
    path_temp = os.path.join(path_output, "temp")

    if not os.path.exists(path_output):
        os.makedirs(path_output)

    if not os.path.exists(path_temp):
        os.makedirs(path_temp)

    # copy input files
    sh.copyfile(file_inv2, os.path.join(path_temp, "inv2.nii"))
    sh.copyfile(file_flash, os.path.join(path_temp, "flash.nii"))

    # convert orig to nifti
    mri_convert(file_orig, os.path.join(path_temp, "orig.nii"))

    # scanner transformation
    scanner_transform(
        os.path.join(path_temp, "inv2.nii"),
        os.path.join(path_temp, "flash.nii"),
        path_temp,
        False,
    )
    scanner_transform(
        os.path.join(path_temp, "flash.nii"),
        os.path.join(path_temp, "inv2.nii"),
        path_temp,
        False,
    )
    scanner_transform(
        os.path.join(path_temp, "inv2.nii"),
        os.path.join(path_temp, "orig.nii"),
        path_temp,
        False,
    )

    # generate coordinate mapping
    generate_coordinate_mapping(
        os.path.join(path_temp, "flash.nii"), 0, path_temp, "flash", False, True
    )

    # scanner transform inv2 to flash
    apply_coordinate_mappings(
        os.path.join(path_temp, "inv2.nii"),  # input
        os.path.join(path_temp, "inv2_2_flash_scanner.nii"),  # cmap
        interpolation="linear",  # nearest or linear
        padding="zero",  # closest, zero or max
        save_data=True,  # save output data to file (boolean)
        overwrite=True,  # overwrite existing results (boolean)
        output_dir=path_temp,  # output directory
        file_name="inv2_apply_scanner",  # base name with file extension for output
    )

    # flirt flash to inv2
    os.chdir(path_temp)
    flirt(
        os.path.join(path_temp, "flash.nii"),
        os.path.join(path_temp, "inv2_apply_scanner_def-img.nii.gz"),
        os.path.join(path_temp, "flash_apply_flirt_def-img.nii.gz"),
        os.path.join(path_temp, "flirt_matrix.txt"),
        cost_func="mutualinfo",
        interp_method="trilinear",
    )

    # apply flirt to flash cmap
    applyxfm = ApplyXFM()
    applyxfm.inputs.in_file = os.path.join(path_temp, "cmap_flash.nii")
    applyxfm.inputs.reference = os.path.join(path_temp, "flash.nii")
    applyxfm.inputs.in_matrix_file = os.path.join(path_temp, "flirt_matrix.txt")
    applyxfm.inputs.interp = "trilinear"
    applyxfm.inputs.padding_size = 0
    applyxfm.inputs.output_type = "NIFTI_GZ"
    applyxfm.inputs.out_file = os.path.join(
        path_temp, "cmap_apply_flirt_def-img.nii.gz"
    )
    applyxfm.inputs.apply_xfm = True
    applyxfm.run()

    # apply scanner transform to flash cmap
    apply_coordinate_mappings(
        os.path.join(path_temp, "cmap_apply_flirt_def-img.nii.gz"),
        os.path.join(path_temp, "flash_2_inv2_scanner.nii"),  # cmap 1
        os.path.join(path_temp, "inv2_2_orig_scanner.nii"),  # cmap 2
        interpolation="linear",  # nearest or linear
        padding="zero",  # closest, zero or max
        save_data=True,  # save output data to file (boolean)
        overwrite=True,  # overwrite existing results (boolean)
        output_dir=path_temp,  # output directory
        file_name="cmap_apply_scanner",  # base name with file extension for output
    )

    # apply deformation to source image
    apply_coordinate_mappings(
        os.path.join(path_temp, "flash.nii"),  # input
        os.path.join(path_temp, "cmap_apply_scanner_def-img.nii.gz"),
        interpolation="linear",  # nearest or linear
        padding="zero",  # closest, zero or max
        save_data=True,  # save output data to file (boolean)
        overwrite=True,  # overwrite existing results (boolean)
        output_dir=path_temp,  # output directory
        file_name="flash_apply_deformation",  # base name with file extension for output
    )

    # rename final deformation examples
    os.rename(
        os.path.join(path_temp, "cmap_apply_scanner_def-img.nii.gz"),
        os.path.join(path_output, "flash2orig.nii.gz"),
    )
    os.rename(
        os.path.join(path_temp, "flash_apply_deformation_def-img.nii.gz"),
        os.path.join(path_output, "flash2orig_example.nii.gz"),
    )

    # clean intermediate files
    if cleanup:
        sh.rmtree(path_temp, ignore_errors=True)
