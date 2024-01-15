# -*- coding: utf-8 -*-
"""Mask brain."""

import os
import shutil as sh

import nibabel as nb
import numpy as np
from nighres.registration import apply_coordinate_mappings, embedded_antsreg
from scipy.ndimage import binary_fill_holes, gaussian_filter
from scipy.ndimage.morphology import binary_dilation

from ..io.filename import get_filename
from ..registration.cmap import expand_coordinate_mapping
from ..registration.transform import scanner_transform

__all__ = ["mask_ana", "clean_ana", "mask_epi"]


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

    Returns
    -------
    None.

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

    Returns
    -------
    None.

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
    ana_reg = apply_coordinate_mappings(
        file_t1,  # input
        file_cmap_reg,  # cmap
        interpolation="linear",  # nearest or linear
        padding="zero",  # closest, zero or max
        save_data=False,
        overwrite=False,
        output_dir=None,
        file_name=None,
    )
    nb.save(ana_reg["result"], file_ana_reg)

    # rigid registration
    embedded_antsreg(
        file_ana_reg,  # source image
        file_epi,  # target image
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
        save_data=True,  # save output data to file
        overwrite=True,  # overwrite existing results
        output_dir=path_t1,  # output directory
        file_name="syn",  # output basename
    )

    # remove unnecessary files
    os.remove(os.path.join(path_t1, "syn_ants-def0.nii.gz"))
    os.remove(os.path.join(path_t1, "syn_ants-invmap.nii.gz"))

    # rename cmap
    os.rename(os.path.join(path_t1, "syn_ants-map.nii.gz"), file_cmap_ants)

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
    cmap_def = apply_coordinate_mappings(
        file_cmap_reg,  # input
        file_cmap_ants,
        interpolation="linear",  # nearest or linear
        padding="zero",  # closest, zero or max
        save_data=False,
        overwrite=False,
        output_dir=None,
        file_name=None,
    )
    nb.save(cmap_def["result"], file_cmap_def)

    # remove outliers and expand
    arr_cmap = cmap_def["result"].get_fdata()

    pts_cmap0 = arr_cmap[0, 0, 0, 0]
    pts_cmap1 = arr_cmap[0, 0, 0, 1]
    pts_cmap2 = arr_cmap[0, 0, 0, 2]

    arr_cmap[arr_cmap == 0] = 0
    arr_cmap[arr_cmap == pts_cmap0] = 0
    arr_cmap[arr_cmap == pts_cmap1] = 0
    arr_cmap[arr_cmap == pts_cmap2] = 0

    output = nb.Nifti1Image(
        arr_cmap, cmap_def["result"].affine, cmap_def["result"].header
    )
    nb.save(output, file_cmap_def)

    expand_coordinate_mapping(
        cmap_in=file_cmap_def,
        path_output=path_t1,
        name_output="cmap_def",
        write_output=True,
    )

    # apply final cmap to t1 and mask
    ana_def = apply_coordinate_mappings(
        file_t1,  # input
        file_cmap_def,
        interpolation="linear",  # nearest or linear
        padding="zero",  # closest, zero or max
        save_data=False,
        overwrite=False,
        output_dir=None,
        file_name=None,
    )
    nb.save(ana_def["result"], file_ana_def)

    mask_def = apply_coordinate_mappings(
        file_mask,  # input
        file_cmap_def,
        interpolation="nearest",  # nearest or linear
        padding="zero",  # closest, zero or max
        save_data=False,
        overwrite=False,
        output_dir=None,
        file_name=None,
    )
    nb.save(mask_def["result"], file_mask_def)

    # finalise mask
    arr_mask = mask_def["result"].get_fdata()
    arr_mask = binary_fill_holes(arr_mask).astype(int)  # fill holes in mask
    arr_mask = binary_dilation(arr_mask, iterations=niter).astype(
        np.float
    )  # dilate mask
    arr_mask = gaussian_filter(arr_mask, sigma=sigma)  # apply gaussian filter

    # write final epi mask
    out_img = nb.Nifti1Image(
        arr_mask, mask_def["result"].affine, mask_def["result"].header
    )
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

    Returns
    -------
    None.

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