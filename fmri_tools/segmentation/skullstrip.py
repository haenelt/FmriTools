# -*- coding: utf-8 -*-
"""Skull stripping tools."""

import os

import nibabel as nb
import numpy as np
from scipy.ndimage.morphology import binary_fill_holes
from scipy.signal import argrelextrema

from ..io.filename import get_filename
from ..registration.transform import apply_header
from ..utils.bias import remove_bias_ants
from ..utils.roi import dilate_fsl, erode_fsl

__all__ = ["skullstrip_flash", "skullstrip_epi", "skullstrip_refined"]


def skullstrip_flash(
    file_in,
    path_output,
    name_output,
    min_val=100,
    max_val=1000,
    flood_fill=False,
    cleanup=False,
):
    """This function computes a brain mask for a partial coverage T2*-weighted
    anatomical image. The mask is used to remove the sagittal sinus during segmentation.
    Thus, the latest echo should be used here to get the largest difference between
    venous and tissue compartments. The brain mask is generated by a simple thresholding
    operation. An intenstiy gradient in posterior-anterior direction is considered by
    thresholding all slices independently. The threshold value is computed from the
    intensity histogram of each slice by getting the minimum of the histogram within a
    predefined range.

    Parameters
    ----------
    file_in : str
        Input path of T2*-weighted anatomy.
    path_output : str
        Path where output is saved.
    name_output : str
        Basename of output file.
    min_val : float, optional
        Minimum threshold of intensity histogram. The default is 100.
    max_val : float, optional
        Maximum threshold of intenstiy histogram. The default is 1000.
    flood_fill : bool, optional
        Apply flood filling of binary mask. The default is False.
    cleanup : bool, optional
        Delete intermediate files. The default is False.

    Returns
    -------
    None.

    """
    # prepare path and filename
    path = os.path.dirname(file_in)
    file = os.path.splitext(os.path.basename(file_in))[0]

    # bias field correction
    remove_bias_ants(
        os.path.join(path, file + ".nii"), os.path.join(path, "b" + file + ".nii")
    )

    # load input data
    mask = nb.load(os.path.join(path, "b" + file + ".nii"))
    mask_array = mask.get_fdata()

    # loop through slices
    for i in range(np.shape(mask_array)[2]):
        # load slice
        temp = np.reshape(mask_array[:, :, i], np.size(mask_array[:, :, i]))

        # get histogram
        bins, edge = np.histogram(temp, 100)
        edge = edge[:-1]

        # find local minimum within defined range
        bin_min = argrelextrema(bins, np.less)
        edge_min = edge[bin_min]

        edge_min[edge_min < min_val] = 0
        edge_min[edge_min > max_val] = 0
        edge_min[edge_min == 0] = np.nan
        edge_min = np.nanmin(edge_min)

        # mask image
        mask_array[:, :, i][mask_array[:, :, i] < edge_min] = 0
        mask_array[:, :, i][mask_array[:, :, i] != 0] = 1

    # flood filling on brain mask
    if flood_fill:
        mask_array = binary_fill_holes(mask_array, structure=np.ones((2, 2, 2)))

    # write output
    output = nb.Nifti1Image(mask_array, mask.affine, mask.header)
    nb.save(output, os.path.join(path_output, name_output + "_mask.nii"))

    # clean intermediate files
    if cleanup:
        os.remove(os.path.join(path, "n4bias.nii"))
        os.remove(os.path.join(path, "b" + file + ".nii"))


def skullstrip_epi(
    file_in, roi_size=5, scale=0.75, nerode=2, ndilate=1, savemask=False, cleanup=True
):
    """Skullstrip input volume by defining an intensity threshold from the inner of the
    brain volume. From a defined mid-point, a brain mask is grown inside the brain. A
    binary filling holes algorithm is applied. To reduce remaining skull within the
    brain mask, the mask is eroded and dilated several times.

    Parameters
    ----------
    file_in : str
        Input file.
    roi_size : int, optional
        Size of cubic roi for image intensity threshold. The default is 5.
    scale : float, optional
        Scale image intensity threshold. The default is 0.75.
    nerode : int, optional
        Number of eroding iterations. The default is 2.
    ndilate : int, optional
        Number of dilating iterations. The default is 1.
    savemask : bool, optional
        Save mask time series. The default is False.
    cleanup : bool, optional
        Delete intermediate files after running. The default is True.

    Returns
    -------
    None.

    """
    # prepare path and filename
    path = os.path.split(file_in)[0]
    file = os.path.split(file_in)[1]

    # load data
    data_img = nb.load(os.path.join(path, file))
    data_array = data_img.get_data()

    # calculate mean intensity
    data_mean = data_array.mean()

    # get point within the brain
    inds = np.transpose(np.nonzero(data_array > data_mean))
    x_mean = np.uint8(np.round((np.max(inds[:, 0]) + np.min(inds[:, 0])) / 2))
    y_mean = np.uint8(np.round((np.max(inds[:, 1]) + np.min(inds[:, 1])) / 2))
    z_mean = np.uint8(np.round((np.max(inds[:, 2]) + np.min(inds[:, 2])) / 2))

    # initialise mask
    mask_array = np.zeros_like(data_array)
    mask_temp_array = np.zeros_like(data_array)
    mask_array[x_mean, y_mean, z_mean] = 1
    mask_temp_array[x_mean, y_mean, z_mean] = 1

    # compute threshold
    roi = data_array[
        np.uint8(np.round(x_mean - roi_size / 2)) : np.uint8(
            np.round(x_mean + roi_size - 1 / 2)
        ),
        np.uint8(np.round(y_mean - roi_size / 2)) : np.uint8(
            np.round(y_mean + roi_size - 1 / 2)
        ),
        np.uint8(np.round(z_mean - roi_size / 2)) : np.uint8(
            np.round(z_mean + roi_size - 1 / 2)
        ),
    ]
    roi_mean = roi.mean()

    # grow mask
    while True:
        coords = np.transpose(np.nonzero(mask_array > 0))

        coords = coords[coords[:, 0] > 0, :]
        coords = coords[coords[:, 1] > 0, :]
        coords = coords[coords[:, 2] > 0, :]
        coords = coords[coords[:, 0] < np.size(data_array, 0) - 1]
        coords = coords[coords[:, 1] < np.size(data_array, 1) - 1]
        coords = coords[coords[:, 2] < np.size(data_array, 2) - 1]

        # calculate neighbour coordinate
        mask_temp_array[coords[:, 0] - 1, coords[:, 1], coords[:, 2]] = 1
        mask_temp_array[coords[:, 0], coords[:, 1] - 1, coords[:, 2]] = 1
        mask_temp_array[coords[:, 0], coords[:, 1], coords[:, 2] - 1] = 1
        mask_temp_array[coords[:, 0] + 1, coords[:, 1], coords[:, 2]] = 1
        mask_temp_array[coords[:, 0], coords[:, 1] + 1, coords[:, 2]] = 1
        mask_temp_array[coords[:, 0], coords[:, 1], coords[:, 2] + 1] = 1

        # delete all old mask elements
        mask_temp_array = mask_temp_array - mask_array
        mask_temp_array[data_array < scale * roi_mean] = 0

        # reinitialise mask_temp
        mask_array = mask_array + mask_temp_array

        # check break condition
        coords = np.transpose(np.nonzero(mask_temp_array == True))
        if len(coords) == 0:
            break

    # flood filling on brain mask
    mask_array = binary_fill_holes(mask_array, structure=np.ones((2, 2, 2)))

    # write mask (intermediate)
    newimg = nb.Nifti1Image(mask_array, data_img.affine, data_img.header)
    newimg.header["dim"][0] = 3
    nb.save(newimg, os.path.join(path, "temp.nii"))

    # erode mask
    for _ in range(nerode):
        erode_fsl(os.path.join(path, "temp.nii"), os.path.join(path, "temp.nii"))

    # dilate mask
    for _ in range(ndilate):
        dilate_fsl(os.path.join(path, "temp.nii"), os.path.join(path, "temp.nii"))

    # load final mask
    temp_img = nb.load(os.path.join(path, "temp.nii"))
    mask_array = temp_img.get_data()

    # write masked image
    data_masked_array = data_array * mask_array
    output = nb.Nifti1Image(data_masked_array, data_img.affine, data_img.header)
    nb.save(output, os.path.join(path, "p" + file))

    # write final output
    if savemask is True:
        newimg = nb.Nifti1Image(mask_array, data_img.affine, data_img.header)
        nb.save(newimg, os.path.join(path, "mask_" + file))

    # clear output
    if cleanup is True:
        os.remove(os.path.join(path, "temp.nii"))


def skullstrip_refined(file_mask1, file_mask2):
    """The purpose of the following function is to enhance the skullstrip mask in native
    space. It uses a second mask which was manually corrected during the freesurfer
    segmentation. This corrected brainmask is converted to native space and multiplied
    with the initial brainmask.

    Parameters
    ----------
    file_mask1 : str
        Brainmask in original space.
    file_mask2 : str
        Manually corrected brainmask in freesurfer space (brain.finalsurfs.mgz).

    Returns
    -------
    file_out : str
        Filename of enhanced brainmask.

    """
    # get output path and basename
    path_output, name_output, _ = get_filename(file_mask1)

    # filename of temporary and enhanced brainmask
    file_temp = os.path.join(path_output, "temp.nii")
    file_out = os.path.join(path_output, name_output + "_enhanced.nii")

    # bring skullstrip_mask from conformed space into original space
    apply_header(
        file_mask2,
        file_mask1,
        file_temp,
        interp_method="nearest",
    )

    # load first brainmask in original space
    mask1 = nb.load(file_mask1)
    mask1_array = mask1.get_fdata()

    # load second brainmask transformed into original space
    mask2 = nb.load(file_temp)
    mask2_array = mask2.get_fdata()

    # make second brainmask binary
    mask2_array[mask2_array == 1] = 0
    mask2_array[mask2_array > 0] = 1

    # multiply both masks
    mask_enhanced_array = mask1_array * mask2_array

    # write enhancec brainmask
    mask_enhanced = nb.Nifti1Image(mask_enhanced_array, mask1.affine, mask1.header)
    nb.save(mask_enhanced, file_out)

    # remove temporary file
    os.remove(file_temp)

    return file_out
