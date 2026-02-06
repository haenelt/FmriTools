# -*- coding: utf-8 -*-
"""Rigid registration."""

import os
import shutil as sh
from pathlib import Path

import nibabel as nb
import numpy as np

from .. import execute_command
from ..io.affine import apply_affine_chunked, read_vox2vox
from ..io.vol import copy_header, mri_convert
from ..segmentation.vol import remove_bias_ants
from ..utils.calc import remove_nans
from .cmap import generate_coordinate_mapping
from .fsl import apply_flirt, flirt
from .transform import apply_coordinate_mapping, scanner_transform

__all__ = ["get_flash2orig", "boundary_based_registration"]


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
    path_output,
    file_mask="",
    wm_bright=False,
    init_reg="header",
    nmax=1000,
    dof=12,
    cleanup=False,
):
    """The script calls the freesurfer bbr method. Inputs are not checked if they are
    valid. Two white surfaces are expected corresponding to left and right hemisphere.
    The target volume is used as freesurfer orig file. The resulting transformation
    matrix is expressed as coordinate mapping and saved as nifti volume. Example files
    are created by applying the resulting coordinate mappings to source and target
    images. The script needs an installation of freesufer.

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
    path_output : str
        Path where output is saved.
    file_mask : str, optional
        File name of mask image in target space. Necessary if init_reg is not "header".
        By default "".
    wm_bright : bool
        WM brighter than GM.
    init_reg : str, optional
        Initialize registration, by default "header".
    nmax : int, optional
        Maximum number of iterations, by default 1000.
    dof : int, optional
        Degrees of freedom (6, 9, or 12), by default 12.
    cleanup : bool, optional
        Delete intermediate files, by default False.
    """
    # set freesurfer path environment
    os.environ["SUBJECTS_DIR"] = path_output

    # freesurfer subject
    sub = "temp"

    # mimic freesurfer folder structure (with some additional folders for intermediate
    # files)
    dir_out = Path(path_output)
    dir_sub = dir_out / sub
    dir_mri = dir_sub / "mri"
    dir_surf = dir_sub / "surf"
    dir_bbr = dir_sub / "bbr"

    # make output folders
    dir_out.mkdir(exist_ok=True, parents=True)
    dir_sub.mkdir(exist_ok=True, parents=True)
    dir_mri.mkdir(exist_ok=True, parents=True)
    dir_surf.mkdir(exist_ok=True, parents=True)
    dir_bbr.mkdir(exist_ok=True, parents=True)

    # change path to output folder
    os.chdir(dir_out)

    # copy surfaces
    sh.copyfile(lh_white, str(dir_surf / "lh.white"))
    sh.copyfile(rh_white, str(dir_surf / "rh.white"))

    # copy volumes
    _file_target = str(dir_mri / "tmp.nii")
    if file_target.endswith(".mgz"):
        mri_convert(file_target, _file_target)
    elif file_target.endswith(".nii"):
        sh.copy(file_target, _file_target)
    else:
        raise ValueError("Invalid file extensions. *.nii or *.mgz file expected!")

    _file_mask = str(dir_mri / "brainmask.mgz")
    if file_mask != "" and file_mask.endswith(".mgz"):
        sh.copy(file_mask, _file_mask)
    elif file_mask != "" and file_mask.endswith(".nii"):
        mri_convert(file_mask, _file_mask)
    else:
        print("No mask detectected. Check inputs if a mask is necessary!")

    ext = "".join(Path(file_source).suffixes)
    sh.copy(file_source, str(dir_bbr / f"source{ext}"))

    # remove nans
    remove_nans(_file_target, str(dir_mri / "orig.nii"))
    remove_nans(str(dir_bbr / f"source{ext}"), str(dir_bbr / f"source{ext}"))

    # scale anatomy
    _orig = nb.load(dir_mri / "orig.nii")
    _arr = _orig.get_fdata()
    q = np.percentile(_arr, [1, 99])
    _arr[_arr < 0.0] = 0.0
    _arr[_arr > q[1]] = q[1]
    _arr /= np.max(_arr)
    _arr *= 4095.0
    output = nb.Nifti1Image(_arr, _orig.affine, _orig.header)
    nb.save(output, dir_mri / "orig.nii")

    # bias field correction
    remove_bias_ants(str(dir_mri / "orig.nii"), str(dir_mri / "borig.nii"))

    # convert orig and brainmask to mgz
    mri_convert(str(dir_mri / "borig.nii"), str(dir_mri / "orig.mgz"))

    # choose initialization method
    if init_reg == "header":
        bbr_var = " --init-header"
    elif init_reg == "freesurfer":
        bbr_var = " --init-coreg"
    elif init_reg == "fsl":
        bbr_var = " --init-fsl"
    else:
        raise ValueError(f"Invalid argument for init_reg: {init_reg}!")

    # choose contrast
    contrast = "--t1" if wm_bright else "--bold"

    # run freesurfer bbr
    os.chdir(dir_bbr)
    command = "bbregister"
    command += f" --s {sub}"
    command += f" --mov {dir_bbr / f'source{ext}'}"
    command += f" {contrast} --reg regheader --gm-proj-abs 1 --wm-proj-abs 1"
    command += f" --nmax {nmax}"
    command += f" --o {str(dir_bbr / 'registered.nii')}"
    command += f" --lta {str(dir_bbr / 'transformation.lta')}"
    command += f" --no-cortex-label --{dof} {bbr_var}"
    command += f" --nocleanup --tmp {dir_bbr}"

    # run
    execute_command(command)

    # get transformation matrix from freesurfer lta file
    M, Minv = read_vox2vox(str(dir_bbr / "transformation.lta"))

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
    nb.save(t2s, dir_out / "target2source.nii.gz")

    # apply cmap to target
    apply_coordinate_mapping(
        _file_target,
        str(dir_out / "target2source.nii.gz"),
        str(dir_out / "target2source_example.nii.gz"),
        interpolation="linear",
    )

    # source to target transformation
    cmap_target = generate_coordinate_mapping(_file_target, pad=0)
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
    s2t_header = copy_header(_file_target)
    s2t_header["dim"][0] = 4
    s2t_header["dim"][4] = 3
    s2t_header["pixdim"][0] = 1
    s2t_header["pixdim"][4] = 1
    s2t_affine = nb.load(_file_target).affine
    s2t = nb.Nifti1Image(arr_cmap_transformed, s2t_affine, s2t_header)
    nb.save(s2t, dir_out / "source2target.nii.gz")

    # apply cmap to source
    apply_coordinate_mapping(
        file_source,
        str(dir_out / "source2target.nii.gz"),
        str(dir_out / "source2target_example.nii.gz"),
        interpolation="linear",
    )

    # clean intermediate files
    if cleanup:
        sh.rmtree(dir_sub, ignore_errors=True)
