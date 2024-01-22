# -*- coding: utf-8 -*-
"""Rigid registration."""

import os
import shutil as sh
import subprocess

import nibabel as nb
import numpy as np

from ..io.affine import apply_affine_chunked, read_vox2vox
from ..io.vol import copy_header, mri_convert
from ..segmentation.mask import clean_ana, mask_ana, mask_epi
from ..utils.bias import remove_bias_ants
from ..utils.calc import multiply_images, remove_nans
from .cmap import generate_coordinate_mapping
from .fsl import apply_flirt, flirt
from .transform import apply_coordinate_mapping, scanner_transform

__all__ = ["get_flash2orig", "boundary_based_registration"]


# parameters for orig skullstrip
NITER_MASK = 3
SIGMA_MASK = 3


def get_flash2orig(file_flash, file_inv2, file_orig, path_output, cleanup=False):
    """This function computes the deformation field for the registration between a
    partial coverage GRE image and the freesurfer orig file. The following steps are
    performed: (1) set output folder structure, (2) get scanner transform GRE <-> inv2
    and inv2 -> orig, (3) generate flash cmap, (4) apply scanner transform inv2 -> GRE,
    (5) get flirt registration GRE -> inv2, (6) apply flirt to GRE cmap, (7) apply
    scanner transform to GRE cmap, (8) apply final deformation to GRE. The function
    needs the FSL environment set.

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
    apply_coordinate_mapping(
        os.path.join(path_temp, "inv2.nii"),
        os.path.join(path_temp, "inv2_2_flash_scanner.nii"),
        os.path.join(path_temp, "inv2_apply_scanner_def-img.nii.gz"),
        interpolation="linear",
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
    apply_flirt(
        os.path.join(path_temp, "cmap_flash.nii"),
        os.path.join(path_temp, "flash.nii"),
        os.path.join(path_temp, "flirt_matrix.txt"),
        os.path.join(path_temp, "cmap_apply_flirt_def-img.nii.gz"),
        "trilinear",
    )

    # concatenate cmaps
    apply_coordinate_mapping(
        os.path.join(path_temp, "flash_2_inv2_scanner.nii"),
        os.path.join(path_temp, "inv2_2_orig_scanner.nii"),
        os.path.join(path_temp, "flash_2_orig_scanner.nii"),
        interpolation="linear",
    )

    # apply scanner transform to flash cmap
    apply_coordinate_mapping(
        os.path.join(path_temp, "cmap_apply_flirt_def-img.nii.gz"),
        os.path.join(path_temp, "flash_2_orig_scanner.nii"),
        os.path.join(path_temp, "cmap_apply_scanner_def-img.nii.gz"),
        interpolation="linear",
    )

    # apply deformation to source image
    apply_coordinate_mapping(
        os.path.join(path_temp, "flash.nii"),  # input
        os.path.join(path_temp, "cmap_apply_scanner_def-img.nii.gz"),
        os.path.join(path_temp, "flash_apply_deformation_def-img.nii.gz"),
        interpolation="linear",  # nearest or linear
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


def boundary_based_registration(
    file_source,
    file_target,
    lh_white,
    rh_white,
    file_ana,
    file_mask,
    file_cmap,
    path_output,
    init_reg="header",
    nmax=1000,
    cleanup=False,
):
    """The script calls the freesurfer bbr method. Inputs are not checked if they are
    valid. Two white surfaces are expected corresponding to left and right
    hemisphere. Therefore, the input should contain a list with two entries. The
    target volume is used as freesurfer orig file. A freesurfer brainmask file is
    generated by skullstripping the target file using a mask from a separate
    anatomical scan. The mask is transformed to orig space by registering ana and
    target. The resulting transformation matrix is expressed as coordinate mapping
    and saved as nifti volume. Example files are created by applying the resulting
    coordinate mappings to source and target images. The script needs an
    installation of freesufer and ants.

    Parameters
    ----------
    file_source : str
        File name of the source image.
    file_target : str
        File name of the target image.
    lh_white : str
        File name of left white surface.
    rh_white : str
        File name of right white surface.
    file_ana : str
        Anatomical image.
    file_mask : str
        Brain mask in anatomical space.
    file_cmap : str
        Coordinate mapping to transform brain mask in anatomica space to the space
        of the source image.
    path_output : str
        Path where output is saved.
    init_reg : str, optional
        Initialize registration, by default "header".
    nmax : int, optional
        Maximum number of iterations, by default 1000.
    cleanup : bool, optional
        Delete intermediate files, by default False.
    """
    # set freesurfer path environment
    os.environ["SUBJECTS_DIR"] = path_output

    # freesurfer subject
    sub = "temp"

    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # mimic freesurfer folder structure (with some additional folders for
    # intermediate files)
    path_sub = os.path.join(path_output, sub)
    path_mri = os.path.join(path_sub, "mri")
    path_surf = os.path.join(path_sub, "surf")
    path_t1 = os.path.join(path_sub, "t1")
    path_bbr = os.path.join(path_sub, "bbr")

    os.makedirs(path_sub)
    os.makedirs(path_mri)
    os.makedirs(path_surf)
    os.makedirs(path_t1)
    os.makedirs(path_bbr)

    # change path to output folder
    os.chdir(path_output)

    # copy surfaces
    sh.copyfile(lh_white, os.path.join(path_surf, "lh.white"))
    sh.copyfile(rh_white, os.path.join(path_surf, "rh.white"))

    # copy volumes
    sh.copy(file_target, os.path.join(path_mri, "orig.nii"))
    sh.copy(file_source, os.path.join(path_bbr, "source.nii"))
    sh.copy(file_ana, os.path.join(path_t1, "T1.nii"))
    sh.copy(file_mask, os.path.join(path_t1, "mask.nii"))

    # remove nans
    remove_nans(os.path.join(path_mri, "orig.nii"), os.path.join(path_mri, "orig.nii"))
    remove_nans(
        os.path.join(path_bbr, "source.nii"), os.path.join(path_bbr, "source.nii")
    )
    remove_nans(os.path.join(path_t1, "T1.nii"), os.path.join(path_t1, "T1.nii"))
    remove_nans(os.path.join(path_t1, "mask.nii"), os.path.join(path_t1, "mask.nii"))

    # get brainmask
    clean_ana(os.path.join(path_t1, "T1.nii"), 1000.0, 4095.0, overwrite=True)
    mask_ana(
        os.path.join(path_t1, "T1.nii"),
        os.path.join(path_t1, "mask.nii"),
        background_bright=False,
    )

    # bias field correction
    remove_bias_ants(
        os.path.join(path_mri, "orig.nii"), os.path.join(path_mri, "borig.nii")
    )

    mask_epi(
        os.path.join(path_mri, "borig.nii"),
        os.path.join(path_t1, "pT1.nii"),
        os.path.join(path_t1, "mask.nii"),
        NITER_MASK,
        SIGMA_MASK,
        file_cmap,
    )

    multiply_images(
        os.path.join(path_mri, "orig.nii"),
        os.path.join(path_t1, "mask_def-img3.nii.gz"),
        os.path.join(path_mri, "brainmask.nii"),
    )

    # convert orig and brainmask to mgz
    mri_convert(os.path.join(path_mri, "orig.nii"), os.path.join(path_mri, "orig.mgz"))
    mri_convert(
        os.path.join(path_mri, "brainmask.nii"), os.path.join(path_mri, "brainmask.mgz")
    )

    # choose initialization method
    if init_reg == "header":
        bbr_var = " --init-header"
    elif init_reg == "freesurfer":
        bbr_var = " --init-coreg"
    elif init_reg == "fsl":
        bbr_var = " --init-fsl"

    # run freesurfer bbr
    os.chdir(path_bbr)
    command = "bbregister"
    command += f" --s {sub}"
    command += f" --mov {os.path.join(path_bbr, 'source.nii')}"
    command += " --bold --reg regheader --gm-proj-abs 1 --wm-proj-abs 1"
    command += f" --nmax {nmax}"
    command += f" --o {os.path.join(path_bbr, 'registered.nii')}"
    command += f" --lta {os.path.join(path_bbr, 'transformation.lta')}"
    command += f" --no-cortex-label --6 {bbr_var}"
    command += f" --nocleanup --tmp {path_bbr}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")

    # get transformation matrix from freesurfer lta file
    M, Minv = read_vox2vox(os.path.join(path_bbr, "transformation.lta"))

    # target to source coordinate mapping
    cmap_source = generate_coordinate_mapping(file_source, pad=0)
    arr_cmap_source = cmap_source.get_fdata()

    xdim = cmap_source.header["dim"][1]
    ydim = cmap_source.header["dim"][2]
    zdim = cmap_source.header["dim"][3]

    x = arr_cmap_source[:, :, :, 0].flatten()
    y = arr_cmap_source[:, :, :, 1].flatten()
    z = arr_cmap_source[:, :, :, 2].flatten()

    source_listed = np.array([x, y, z]).T
    source_transformed = apply_affine_chunked(M, source_listed)

    x_new = np.reshape(source_transformed[:, 0], (xdim, ydim, zdim))
    y_new = np.reshape(source_transformed[:, 1], (xdim, ydim, zdim))
    z_new = np.reshape(source_transformed[:, 2], (xdim, ydim, zdim))

    arr_cmap_transformed = np.zeros_like(arr_cmap_source)
    arr_cmap_transformed[:, :, :, 0] = x_new
    arr_cmap_transformed[:, :, :, 1] = y_new
    arr_cmap_transformed[:, :, :, 2] = z_new

    # nibabel instance of final cmap
    t2s_header = copy_header(file_source)
    t2s_header["dim"][0] = 4
    t2s_header["dim"][4] = 3
    t2s_header["pixdim"][0] = 1
    t2s_header["pixdim"][4] = 1
    t2s_affine = nb.load(file_source).affine
    t2s = nb.Nifti1Image(arr_cmap_transformed, t2s_affine, t2s_header)

    # apply cmap to target
    apply_coordinate_mapping(
        file_target,
        t2s.get_fdata(),
        os.path.join(path_output, "target2source_example.nii.gz"),
        interpolation="linear",
    )

    # source to target transformation
    cmap_target = generate_coordinate_mapping(file_target, pad=0)
    arr_cmap_target = cmap_target.get_fdata()

    xdim = cmap_target.header["dim"][1]
    ydim = cmap_target.header["dim"][2]
    zdim = cmap_target.header["dim"][3]

    # transform source volume
    x = arr_cmap_target[:, :, :, 0].flatten()
    y = arr_cmap_target[:, :, :, 1].flatten()
    z = arr_cmap_target[:, :, :, 2].flatten()

    target_listed = np.array([x, y, z]).T
    target_transformed = apply_affine_chunked(Minv, target_listed)

    x_new = np.reshape(target_transformed[:, 0], (xdim, ydim, zdim))
    y_new = np.reshape(target_transformed[:, 1], (xdim, ydim, zdim))
    z_new = np.reshape(target_transformed[:, 2], (xdim, ydim, zdim))

    arr_cmap_transformed = np.zeros_like(arr_cmap_target)
    arr_cmap_transformed[:, :, :, 0] = x_new
    arr_cmap_transformed[:, :, :, 1] = y_new
    arr_cmap_transformed[:, :, :, 2] = z_new

    # nibabel instance of final cmap
    s2t_header = copy_header(file_target)
    s2t_header["dim"][0] = 4
    s2t_header["dim"][4] = 3
    s2t_header["pixdim"][0] = 1
    s2t_header["pixdim"][4] = 1
    s2t_affine = nb.load(file_target).affine
    s2t = nb.Nifti1Image(arr_cmap_transformed, s2t_affine, s2t_header)

    # apply cmap to source
    apply_coordinate_mapping(
        file_source,
        s2t.get_fdata(),
        os.path.join(path_output, "source2target_example.nii.gz"),
        interpolation="linear",
    )

    # write output
    nb.save(t2s, os.path.join(path_output, "target2source.nii.gz"))
    nb.save(s2t, os.path.join(path_output, "source2target.nii.gz"))

    # clean intermediate files
    if cleanup:
        sh.rmtree(path_sub, ignore_errors=True)
