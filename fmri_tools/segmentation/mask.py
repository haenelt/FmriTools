# -*- coding: utf-8 -*-
"""Mask brain."""

import os
import shutil as sh
import subprocess

import nibabel as nb
import numpy as np
from scipy.ndimage import (binary_dilation, binary_erosion, binary_fill_holes,
                           gaussian_filter)

from ..io.filename import get_filename
from ..io.vol import mri_convert
from ..matlab import MatlabCommand
from ..registration.cmap import expand_coordinate_mapping
from ..registration.nonrigid import embedded_antsreg
from ..registration.transform import (apply_coordinate_mapping,
                                      scanner_transform)
from ..segmentation.skullstrip import skullstrip_bet

__all__ = [
    "mask_ana",
    "clean_ana",
    "mask_epi",
    "deweight_mask",
    "mask_nuisance",
    "mask_ribbon",
]


def mask_ana(t1, mask, background_bright=False):
    """This function masked an image with a corresponding binary mask by multiplication.
    The masked volume is saved in the same folder as the input image with the prefix p.

    Parameters
    ----------
    t1 : str
        Input anatomy.
    mask : str
        Corresponding binary mask.
    background_bright : bool, optional
        Set values outside mask to maximum value. The default is False.
    """
    # get path and filename of anatomy
    path_t1 = os.path.dirname(t1)
    if os.path.splitext(os.path.basename(t1))[1] == ".gz":
        name_t1 = os.path.splitext(os.path.splitext(os.path.basename(t1))[0])[0]
    else:
        name_t1 = os.path.splitext(os.path.basename(t1))[0]

    # load anatomy
    ana_img = nb.load(t1)
    ana_array = ana_img.get_fdata()

    # load mask
    mask_img = nb.load(mask)
    mask_array = mask_img.get_fdata()

    # multiply images
    masked_ana_array = ana_array * mask_array

    # set all outside mask values to maximum to mimic bright CSF values
    if background_bright:
        masked_ana_array[mask_array == 0] = np.max(ana_array)

    # write masked anatomy
    out_img = nb.Nifti1Image(masked_ana_array, ana_img.affine, ana_img.header)
    nb.save(out_img, os.path.join(path_t1, "p" + name_t1 + ".nii"))


def mask_epi(file_epi, file_t1, file_mask, niter, sigma, file_reg=""):
    """This function masks a mean epi image based on a skullstrip mask of the
    corresponding anatomy. The mask is transformed to native epi space via an initial
    transformation or via scanner coordinates. A rigid registration is applied to ensure
    a match between mask and epi. Finally, holes in the mask are filled, the mask is
    dilated and a Gaussian filter is applied. The masked epi is saved in the same folder
    with the prefix p.

    Parameters
    ----------
    file_epi : str
        Input mean epi image.
    file_t1 : str
        Input of corresponding skullstripped anatomy.
    file_mask : str
        Input of skullstrip mask of the corresponding anatomy.
    niter : int
        Number of dilation iterations.
    sigma : float
        Gaussian smoothing kernel.
    file_reg : str, optional
        Filename of ana -> epi coordinate mapping. The default is "".
    """
    # get paths and filenames
    path_t1, name_t1, _ = get_filename(file_t1)
    path_epi, name_epi, _ = get_filename(file_epi)

    if file_reg:
        _, _, ext_reg = get_filename(file_reg)
    else:
        ext_reg = ".nii.gz"

    # filenames
    file_cmap_reg = os.path.join(path_t1, "cmap_reg" + ext_reg)
    file_cmap_ants = os.path.join(path_t1, "cmap_ants.nii.gz")
    file_cmap_def = os.path.join(path_t1, "cmap_def.nii.gz")
    file_ana_reg = os.path.join(path_t1, "ana_reg.nii.gz")
    file_ana_def = os.path.join(path_t1, "ana_def.nii.gz")
    file_mask_def = os.path.join(path_t1, "mask_def.nii.gz")
    file_mask_def2 = os.path.join(path_t1, "mask_def2.nii.gz")

    # get initial ana -> epi transformation from existing cmap or header
    if file_reg:
        sh.copyfile(file_reg, file_cmap_reg)
    else:
        scanner_transform(file_t1, file_epi, path_t1, True)
        os.rename(
            os.path.join(path_t1, name_t1 + "_2_" + name_epi + "_scanner.nii.gz"),
            file_cmap_reg,
        )

    # scanner transform peeled t1 to epi
    apply_coordinate_mapping(
        file_t1, file_cmap_reg, file_ana_reg, interpolation="linear"
    )

    # rigid registration
    embedded_antsreg(
        file_ana_reg,  # source image
        file_epi,  # target image
        output_dir=path_t1,  # output directory
        run_rigid=True,  # whether or not to run a rigid registration first
        rigid_iterations=1000,  # number of iterations in the rigid step
        run_affine=False,  # whether or not to run an affine registration first
        affine_iterations=0,  # number of iterations in the affine step
        run_syn=False,  # whether or not to run a SyN registration
        coarse_iterations=0,  # number of iterations at the coarse level
        medium_iterations=0,  # number of iterations at the medium level
        fine_iterations=0,  # number of iterations at the fine level
        cost_function="CrossCorrelation",  # CrossCorrelation or MutualInformation
        interpolation="Linear",  # interpolation for registration result (NeareastNeighbor or Linear)
        convergence=1e-6,  # threshold for convergence (can make algorithm very slow)
        ignore_affine=True,  # ignore the affine matrix information extracted from the image header
        ignore_header=True,  # ignore the orientation information and affine matrix information extracted from the image header
    )

    # remove unnecessary files
    _, name_ana_reg, _ = get_filename(file_ana_reg)
    os.remove(os.path.join(path_t1, f"{name_ana_reg}_ants-def.nii.gz"))
    os.remove(os.path.join(path_t1, f"{name_ana_reg}_ants-invmap.nii.gz"))

    # rename cmap
    os.rename(os.path.join(path_t1, f"{name_ana_reg}_ants-map.nii.gz"), file_cmap_ants)

    # remove outliers and expand
    cmap = nb.load(file_cmap_ants)
    arr_cmap = cmap.get_fdata()

    pts_cmap0 = arr_cmap[0, 0, 0, 0]
    pts_cmap1 = arr_cmap[0, 0, 0, 1]
    pts_cmap2 = arr_cmap[0, 0, 0, 2]

    arr_cmap[arr_cmap == 0] = 0
    arr_cmap[arr_cmap == pts_cmap0] = 0
    arr_cmap[arr_cmap == pts_cmap1] = 0
    arr_cmap[arr_cmap == pts_cmap2] = 0

    output = nb.Nifti1Image(arr_cmap, cmap.affine, cmap.header)
    nb.save(output, file_cmap_ants)

    expand_coordinate_mapping(
        cmap_in=file_cmap_ants,
        path_output=path_t1,
        name_output="cmap_ants",
        write_output=True,
    )

    # apply ants cmap to header transformation
    cmap_def = apply_coordinate_mapping(
        file_cmap_reg, file_cmap_ants, file_cmap_def, interpolation="linear"
    )

    # remove outliers and expand
    arr_cmap = cmap_def.get_fdata()

    pts_cmap0 = arr_cmap[0, 0, 0, 0]
    pts_cmap1 = arr_cmap[0, 0, 0, 1]
    pts_cmap2 = arr_cmap[0, 0, 0, 2]

    arr_cmap[arr_cmap == 0] = 0
    arr_cmap[arr_cmap == pts_cmap0] = 0
    arr_cmap[arr_cmap == pts_cmap1] = 0
    arr_cmap[arr_cmap == pts_cmap2] = 0

    output = nb.Nifti1Image(arr_cmap, cmap_def.affine, cmap_def.header)
    nb.save(output, file_cmap_def)

    expand_coordinate_mapping(
        cmap_in=file_cmap_def,
        path_output=path_t1,
        name_output="cmap_def",
        write_output=True,
    )

    # apply final cmap to t1 and mask
    apply_coordinate_mapping(
        file_t1, file_cmap_def, file_ana_def, interpolation="linear"
    )
    mask_def = apply_coordinate_mapping(
        file_mask, file_cmap_def, file_mask_def, interpolation="nearest"
    )

    # finalise mask
    arr_mask = mask_def.get_fdata()
    arr_mask = binary_fill_holes(arr_mask).astype(int)  # fill holes in mask
    arr_mask = binary_dilation(arr_mask, iterations=niter).astype(
        np.float64
    )  # dilate mask
    arr_mask = gaussian_filter(arr_mask, sigma=sigma)  # apply gaussian filter

    # write final epi mask
    out_img = nb.Nifti1Image(arr_mask, mask_def.affine, mask_def.header)
    nb.save(out_img, file_mask_def2)

    # multiply epi and binary mask
    epi_img = nb.load(file_epi)
    arr_epi = epi_img.get_fdata()
    arr_epi *= arr_mask  # multiply epi and mask

    # write masked epi
    out_img = nb.Nifti1Image(arr_epi, epi_img.affine, epi_img.header)
    nb.save(out_img, os.path.join(path_epi, "p" + name_epi + ".nii"))


def clean_ana(file_in, min_value, new_range, overwrite=True):
    """This function removes ceiling values from a computed T1 map of an mp2rage
    acquisition. Low intensity values are removed and the data range is normalised to a
    defined new range. The input file should be either a nifti or a compressed nifti
    file.

    Parameters
    ----------
    file_in : str
        Filename of input image.
    min_value : float
        Threshold of low intensity values.
    new_range : float
        Arbitrary new data range.
    overwrite : bool, optional
        If set, the input image is overwritten with the cleaned data set. The
        default is True.
    """
    # load data
    data = nb.load(file_in)
    data_array = data.get_fdata()

    # remove ceiling
    data_array[data_array == 0] = np.max(data_array)

    # remove low intensity values
    data_array[data_array <= min_value] = 0
    data_array = data_array - np.min(data_array[data_array != 0])
    data_array[data_array <= 0] = 0

    # normalise to new data range
    data_array = data_array / np.max(data_array) * new_range

    # output cleaned dataset
    output = nb.Nifti1Image(data_array, data.affine, data.header)
    if overwrite:
        nb.save(output, file_in)
    else:
        path = os.path.dirname(file_in)
        if os.path.splitext(file_in)[1] == ".gz":
            basename = os.path.splitext(os.path.splitext(os.path.basename(file_in))[0])[
                0
            ]
            nb.save(output, os.path.join(path, basename + "_clean.nii.gz"))
        else:
            basename = os.path.splitext(os.path.basename(file_in))[0]
            nb.save(output, os.path.join(path, basename + "_clean.nii"))


def deweight_mask(
    file_in,
    mask_in,
    mask_max=0.25,
    sigma_gaussian=10.0,
    write_output=False,
    path_output=None,
):
    """This function computes a binary mask by pooling all voxels above a given
    threshold and replaces all image voxels by its gaussian filtered image voxels within
    this binary mask.

    Parameters
    ----------
    file_in : str
        Filename of input image.
    mask_in : str
        Filename of input mask.
    mask_max : float, optional
        Cutoff threshold. The default is 0.25.
    sigma_gaussian : float, optional
        Sigma for gaussian filter. The default is 10.0.
    write_output : bool, optional
        Write output image The default is None.
    path_output : str, optional
        Path where output is written. The default is None.

    Returns
    -------
    data_array : np.ndarray
        Image matrix with filtered voxels.

    """
    # get basename of phase file
    _, name_file, ext_file = get_filename(file_in)

    # load unwrapped phase data
    data = nb.load(file_in)
    data_array = data.get_fdata()

    # load standard deviation data
    mask = nb.load(mask_in)
    mask_array = mask.get_fdata()

    # threshold standard deviation
    mask_array[mask_array < mask_max] = 0
    mask_array[mask_array != 0] = 1

    # apply gaussian filter to phase data
    data_array_gaussian = gaussian_filter(
        data_array,
        sigma_gaussian,
        order=0,
        output=None,
        mode="reflect",
        cval=0.0,
        truncate=4.0,
    )

    # replace data
    data_array[mask_array == 1] = data_array_gaussian[mask_array == 1]

    # write output
    if write_output:
        output = nb.Nifti1Image(data_array, data.affine, data.header)
        nb.save(output, os.path.join(path_output, name_file + "_filtered" + ext_file))

    return data_array


def mask_nuisance(
    orig_in,
    deformation,
    path_output,
    nerode_white=1,
    nerode_csf=1,
    segmentation=True,
    cleanup=True,
):
    """This function calculates WM and CSF masks in space of the functional time series.
    It uses SPM to compute WM and CSF probability maps. These maps are masked with a
    skullstrip mask and transformed to native epi space.

    Parameters
    ----------
    orig_in : str
        Input anatomy (orig.mgz).
    deformation : str
        Coordinate mapping for ana to epi transformation.
    path_output : str
        Path where output is saved.
    nerode_white : int, optional
        Number of wm mask eroding steps. The default is 1.
    nerode_csf : int, optional
        Number of csf mask eroding steps. The default is 1.
    segmentation : bool, optional
        Do not calculate new masks to not rerun everything. The default is True.
    cleanup : bool, optional
        Delete intermediate files. The default is True.
    """
    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # get filename without file extension of input file
    file = os.path.splitext(os.path.basename(orig_in))[0]

    # convert to nifti format
    mri_convert(orig_in, os.path.join(path_output, file + ".nii"))

    # bet skullstrip mask
    skullstrip_bet(
        os.path.join(path_output, file + ".nii"), os.path.join(path_output, "_mask.nii")
    )

    # segmentation
    if segmentation:
        matlab = MatlabCommand(
            "ft_skullstrip_spm12", os.path.join(path_output, file + ".nii"), path_output
        )
        matlab.run()

    # load tissue maps
    wm_array = nb.load(
        os.path.join(path_output, "skull", "c2" + file + ".nii")
    ).get_fdata()
    csf_array = nb.load(
        os.path.join(path_output, "skull", "c3" + file + ".nii")
    ).get_fdata()
    mask_array = nb.load(os.path.join(path_output, "bet_mask.nii")).get_fdata()

    # binarize
    wm_array[wm_array > 0] = 1
    csf_array[csf_array > 0] = 1

    # apply brain mask
    wm_array = wm_array * mask_array
    csf_array = csf_array * mask_array

    # erode wm
    wm_array = binary_erosion(
        wm_array,
        structure=None,
        iterations=nerode_white,
        mask=None,
        output=None,
        border_value=0,
        origin=0,
        brute_force=False,
    )

    # erode csf
    csf_array = binary_erosion(
        csf_array,
        structure=None,
        iterations=nerode_csf,
        mask=None,
        output=None,
        border_value=0,
        origin=0,
        brute_force=False,
    )

    # write wm and csf mask
    data_img = nb.load(orig_in)
    wm_out = nb.Nifti1Image(wm_array, data_img.affine, data_img.header)
    nb.save(wm_out, os.path.join(path_output, "wm_mask_orig.nii"))
    csf_out = nb.Nifti1Image(csf_array, data_img.affine, data_img.header)
    nb.save(csf_out, os.path.join(path_output, "csf_mask_orig.nii"))

    # apply deformation to mask
    apply_coordinate_mapping(
        os.path.join(path_output, "wm_mask_orig.nii"),
        deformation,
        os.path.join(path_output, "wm_mask.nii.gz"),
        interpolation="nearest",
    )

    apply_coordinate_mapping(
        os.path.join(path_output, "csf_mask_orig.nii"),
        deformation,
        os.path.join(path_output, "csf_mask.nii.gz"),
        interpolation="nearest",
    )

    # cleanup
    if cleanup:
        os.remove(os.path.join(path_output, "bet_mask.nii"))
        os.remove(os.path.join(path_output, "csf_mask_orig.nii"))
        os.remove(os.path.join(path_output, "wm_mask_orig.nii"))
        os.remove(os.path.join(path_output, "orig.nii"))
        sh.rmtree(os.path.join(path_output, "skull"), ignore_errors=True)


def mask_ribbon(subjects_dir, sub):
    """Computes a ribbon mask, based on a freesurfer segmentation.

    Parameters
    ----------
    subjects_dir : str
        Path of subjects directory.
    sub : str
        Subjects name.
    """
    left_whitelabel = 2
    left_ribbonlabel = 3
    right_whitelabel = 41
    right_ribbonlabel = 42

    # set freesurfer path environment
    os.environ["SUBJECTS_DIR"] = subjects_dir

    command = "mris_volmask"
    command += f" --label_left_ribbon {left_ribbonlabel}"
    command += f" --label_left_white {left_whitelabel}"
    command += f" --label_right_ribbon {right_ribbonlabel}"
    command += f" --label_right_white {right_whitelabel}"
    command += f" --save_ribbon {sub}"

    print("Execute: " + command)
    try:
        subprocess.run([command], shell=True, check=False)
    except subprocess.CalledProcessError:
        print("Execuation failed!")
