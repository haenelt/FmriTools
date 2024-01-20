# -*- coding: utf-8 -*-
"""Coordinate mapping representing image deformations."""

import os
import random

import nibabel as nb
import numpy as np
from numpy.matlib import repmat

from ..io.affine import apply_affine_chunked
from ..io.filename import get_filename

__all__ = [
    "clean_coordinate_mapping",
    "crop_coordinate_mapping",
    "generate_coordinate_mapping",
    "expand_coordinate_mapping",
    "remove_edge_cmap",
]


def clean_coordinate_mapping(
    cmap_source, cmap_target, overwrite_file=True, save_mask=False
):
    """Voxels in the target coordinate mapping are masked out based on found voxel
    displacements in the source coordinate mapping. This is done to remove smeared
    regions caused by interpolations with background values in the case of deforming a
    slab within a larger image array.

    Parameters
    ----------
    cmap_source : str
        Filename of source coordinate mapping.
    cmap_target : str
        Filename of target coordinate mapping.
    overwrite_file : bool, optional
        Overwrite target coordinate mapping. The default is True.
    save_mask : bool, optional
        Write out mask. The default is False.

    Returns
    -------
    results : dict
        Cleaned cmap and mask.

    """
    # get filename
    path_file, _, _ = get_filename(cmap_target)

    # load data
    cmap1_img = nb.load(cmap_source)
    cmap1_array = cmap1_img.get_fdata()

    cmap2_img = nb.load(cmap_target)
    cmap2_array = cmap2_img.get_fdata()

    mask_img = nb.load(cmap_target)
    mask_img.header["dim"][0] = 3
    mask_img.header["dim"][4] = 1
    mask_array = np.zeros_like(mask_img.get_fdata()[:, :, :, 0])

    x_max = cmap2_img.header["dim"][1]
    y_max = cmap2_img.header["dim"][2]
    z_max = cmap2_img.header["dim"][3]

    # get nearest voxel coordinates
    x0 = np.floor(cmap1_array[:, :, :, 0].flatten()).astype(int)
    x1 = np.ceil(cmap1_array[:, :, :, 0].flatten()).astype(int)
    y0 = np.floor(cmap1_array[:, :, :, 1].flatten()).astype(int)
    y1 = np.ceil(cmap1_array[:, :, :, 1].flatten()).astype(int)
    z0 = np.floor(cmap1_array[:, :, :, 2].flatten()).astype(int)
    z1 = np.ceil(cmap1_array[:, :, :, 2].flatten()).astype(int)

    # exclude voxels which do not fit in the target array
    outlier = []
    outlier.extend(np.where(x0 < 0)[0])
    outlier.extend(np.where(x1 < 0)[0])
    outlier.extend(np.where(y0 < 0)[0])
    outlier.extend(np.where(y1 < 0)[0])
    outlier.extend(np.where(z0 < 0)[0])
    outlier.extend(np.where(z1 < 0)[0])
    outlier.extend(np.where(x0 >= x_max)[0])
    outlier.extend(np.where(x1 >= x_max)[0])
    outlier.extend(np.where(y0 >= y_max)[0])
    outlier.extend(np.where(y1 >= y_max)[0])
    outlier.extend(np.where(z0 >= z_max)[0])
    outlier.extend(np.where(z1 >= z_max)[0])
    outlier = list(np.unique(outlier))

    x0 = np.delete(x0, outlier)
    x1 = np.delete(x1, outlier)
    y0 = np.delete(y0, outlier)
    y1 = np.delete(y1, outlier)
    z0 = np.delete(z0, outlier)
    z1 = np.delete(z1, outlier)

    # get final mask
    mask_array[x0, y0, z0] = 1
    mask_array[x1, y1, z1] = 1
    mask_array[x1, y0, z0] = 1
    mask_array[x0, y1, z0] = 1
    mask_array[x0, y0, z1] = 1
    mask_array[x1, y1, z0] = 1
    mask_array[x0, y1, z1] = 1
    mask_array[x1, y0, z1] = 1

    # apply mask to cmap
    cmap2_array[:, :, :, 0] *= mask_array
    cmap2_array[:, :, :, 1] *= mask_array
    cmap2_array[:, :, :, 2] *= mask_array

    # get output
    results = dict()
    results["cmap"] = nb.Nifti1Image(cmap2_array, cmap2_img.affine, cmap2_img.header)
    results["mask"] = nb.Nifti1Image(mask_array, mask_img.affine, mask_img.header)

    # write output
    if overwrite_file:
        nb.save(results["cmap"], cmap_target)

    if save_mask:
        nb.save(results["mask"], os.path.join(path_file, "cmap_mask.nii"))

    return results


def crop_coordinate_mapping(file_in, pad=0, overwrite_file=True, path_output=""):
    """Crops a padded coordinate mapping. The output file can either overwrite the input
    file or a new file is created with a suffix in a defined output directory.

    Parameters
    ----------
    file_in : str
        File name of input file.
    pad : int, optional
        Image padding size.. The default is 0.
    overwrite_file : bool, optional
        Output file overwrites input file.. The default is True.
    path_output : str, optional
        Path where output is saved if input file is not overwritten. The
        default is "".

    Returns
    -------
    None.

    """
    # define output folder
    if path_output is not None:
        if not os.path.exists(path_output):
            os.makedirs(path_output)

    # get input path and file name
    _, file, ext = get_filename(file_in)

    # load data
    data_img = nb.load(file_in)
    data_array = data_img.get_fdata()

    # get matrix size
    x_size = np.size(data_array, 0)
    y_size = np.size(data_array, 1)
    z_size = np.size(data_array, 2)

    # crop image matrix
    data_array = data_array[
        pad : x_size - pad, pad : y_size - pad, pad : z_size - pad, :
    ]

    # write cropped coordinate mapping
    output = nb.Nifti1Image(data_array, data_img.affine, data_img.header)
    output.set_data_dtype(np.float)

    # write coordinate mapping for each time point
    if overwrite_file is True:
        os.remove(file_in)
        nb.save(output, file_in)
    else:
        file_out = os.path.join(path_output, file + "_crop" + ext)
        nb.save(output, file_out)


def generate_coordinate_mapping(
    file_in, pad, path_output=None, suffix=None, time=False, write_output=False
):
    """Generates coordinate mapping for an input volume. Either one or multiple
    coordinate maps are saved in the output folder depending on the dimensionality
    (3d or 4d) of the input image. Image padding can be applied which expands the image
    matrix of each axis in both directions.

    Parameters
    ----------
    file_in : str
        File name of input file.
    pad : int
        Image padding size.
    path_output : str, optional
        Path where output is saved. The default is None.
    suffix : str, optional
        Add suffix to file name. The default is None.
    time : bool, optional
        Compute coordinate map for each time step. The default is False.
    write_output : bool, optional
        Write nifti volume. The default is False.

    Returns
    -------
    output : niimg
        Coordinate mapping.

    """
    # create output folder
    if path_output:
        if not os.path.exists(path_output):
            os.makedirs(path_output)

    # load data
    data_img = nb.load(file_in)

    # get matrix size
    x_size = data_img.header["dim"][1] + 2 * pad
    y_size = data_img.header["dim"][2] + 2 * pad
    z_size = data_img.header["dim"][3] + 2 * pad

    if time is False:
        t_size = 1
    else:
        t_size = data_img.header["dim"][4]

    # define coordinate
    coordinate_mapping = np.zeros((x_size, y_size, z_size, 3), dtype="float")

    # coordinate mapping in x-direction
    x = np.array(np.arange(-pad, x_size - pad, 1), dtype="float")
    x = np.transpose(repmat(x, y_size, 1))
    x = np.dstack([x] * z_size)

    # coordinate mapping in y-direction
    y = np.array(np.arange(-pad, y_size - pad), dtype="float")
    y = repmat(y, x_size, 1)
    y = np.dstack([y] * z_size)

    # coordinate mapping in z-direction
    z = np.ones((x_size, y_size, z_size))
    z = np.arange(-pad, z_size - pad) * z

    # merge directions
    coordinate_mapping[:, :, :, 0] = x
    coordinate_mapping[:, :, :, 1] = y
    coordinate_mapping[:, :, :, 2] = z

    # write coordinate mapping
    output = nb.Nifti1Image(coordinate_mapping, data_img.affine, nb.Nifti1Header())

    # write coordinate mapping for each time point
    if t_size == 1:
        if write_output:
            file_out = os.path.join(path_output, "cmap_" + suffix + ".nii")
            nb.save(output, file_out)
    else:
        for i in range(t_size):
            if write_output:
                file_out = os.path.join(
                    path_output, "cmap_" + suffix + "_" + str(i) + ".nii"
                )
                nb.save(output, file_out)

    return output


def expand_coordinate_mapping(
    cmap_in, path_output=None, name_output=None, write_output=False
):
    """This function removes black background in a coordinate mapping to omit
    interpolation problems at the edges of a coordinate slab within a larger volume.
    Based on the cmap, a transformation matrix is computed from randomly sampled data
    points within the slab. The transformation matrix is then applied to all background
    voxels. Hence, this method is only really precise for coordinate mappings
    representing an affine transformation. However, this function can also be applied to
    nonlinear coordinate mappings since the preliminary goal is to avoid problems at the
    slab edges. Therefore, the actual data sampling should not be affected. The code
    snippet for computing the transformation matrix is taken from [1].

    Parameters
    ----------
    cmap_in : str
        Filename of coordinate mapping.
    path_output : str, optional
        Path where output is written. The default is None.
    name_output : str, optional
        Basename of output volume. The default is None.
    write_output : bool, optional
        Write nifti volume. The default is False.

    Returns
    -------
    output : niimg
        Corrected coordinate mapping.

    References
    -------
    .. [1] https://stackoverflow.com/questions/56220626/how-to-compute-
    conformal-affine-transformation

    """
    # get file extension of cmap
    _, _, ext_cmap = get_filename(cmap_in)

    # load target cmap and generate source cmap
    cmap_target = nb.load(cmap_in)
    cmap_source = generate_coordinate_mapping(cmap_in, pad=0)

    arr_cmap_target = cmap_target.get_fdata()
    arr_cmap_source = cmap_source.get_fdata()

    # get image dimensions
    xdim = cmap_source.header["dim"][1]
    ydim = cmap_source.header["dim"][2]
    zdim = cmap_source.header["dim"][3]

    # random selection of 4 data points
    s_coords = []
    t_coords = []

    pts = np.where(arr_cmap_target[:, :, :, 0] != 0)
    r = random.sample(np.arange(len(pts[0])).tolist(), len(pts[0]))[:4]

    s_coords.append([pts[0][r[0]], pts[1][r[0]], pts[2][r[0]]])
    s_coords.append([pts[0][r[1]], pts[1][r[1]], pts[2][r[1]]])
    s_coords.append([pts[0][r[2]], pts[1][r[2]], pts[2][r[2]]])
    s_coords.append([pts[0][r[3]], pts[1][r[3]], pts[2][r[3]]])

    t_coords.append(
        arr_cmap_target[s_coords[0][0], s_coords[0][1], s_coords[0][2], :].tolist()
    )
    t_coords.append(
        arr_cmap_target[s_coords[1][0], s_coords[1][1], s_coords[1][2], :].tolist()
    )
    t_coords.append(
        arr_cmap_target[s_coords[2][0], s_coords[2][1], s_coords[2][2], :].tolist()
    )
    t_coords.append(
        arr_cmap_target[s_coords[3][0], s_coords[3][1], s_coords[3][2], :].tolist()
    )

    # get transformation matrix (target -> source)
    t_length = len(t_coords)
    B = np.vstack([np.transpose(t_coords), np.ones(t_length)])
    D = 1.0 / np.linalg.det(B)

    M = [
        [
            (-1) ** i * D * np.linalg.det(np.delete(np.vstack([R, B]), (i + 1), axis=0))
            for i in range(t_length)
        ]
        for R in np.transpose(s_coords)
    ]
    arr_rot, vec_trans = np.hsplit(np.array(M), [t_length - 1])
    vec_trans = np.transpose(vec_trans)[0]

    # unittests
    print("Test cmap expansion:")
    for pt, ps in zip(np.array(t_coords), np.array(s_coords)):
        image_p = np.dot(arr_rot, pt) + vec_trans
        result = "[OK]" if np.allclose(image_p, ps) else "[ERROR]"
        print(pt, " mapped to: ", image_p, " ; expected: ", ps, result)

    # get affine transformation matrix by adding translation vector
    trans_forward = np.zeros((4, 4))
    trans_forward[:3, :3] = arr_rot
    trans_forward[:3, -1] = vec_trans
    trans_forward[-1, -1] = 1

    # get final transformation matrix (source -> target)
    trans_inverse = np.linalg.inv(trans_forward)

    # transform source volume
    x = arr_cmap_source[:, :, :, 0].flatten()
    y = arr_cmap_source[:, :, :, 1].flatten()
    z = arr_cmap_source[:, :, :, 2].flatten()

    source_listed = np.array([x, y, z]).T
    source_transformed = apply_affine_chunked(trans_inverse, source_listed)

    x_new = np.reshape(source_transformed[:, 0], (xdim, ydim, zdim))
    y_new = np.reshape(source_transformed[:, 1], (xdim, ydim, zdim))
    z_new = np.reshape(source_transformed[:, 2], (xdim, ydim, zdim))

    # overwrite new cmap with old coordinate mapping (so that only background remains)
    x_new[arr_cmap_target[:, :, :, 0] > 0] = arr_cmap_target[:, :, :, 0][
        arr_cmap_target[:, :, :, 0] > 0
    ]
    y_new[arr_cmap_target[:, :, :, 1] > 0] = arr_cmap_target[:, :, :, 1][
        arr_cmap_target[:, :, :, 1] > 0
    ]
    z_new[arr_cmap_target[:, :, :, 2] > 0] = arr_cmap_target[:, :, :, 2][
        arr_cmap_target[:, :, :, 2] > 0
    ]

    # overwrite input cmap array with final cmap array
    arr_cmap_target[:, :, :, 0] = x_new
    arr_cmap_target[:, :, :, 1] = y_new
    arr_cmap_target[:, :, :, 2] = z_new

    # nibabel instance of final cmap
    output = nb.Nifti1Image(arr_cmap_target, cmap_target.affine, cmap_target.header)

    # write output
    if write_output:
        nb.save(output, os.path.join(path_output, name_output + ext_cmap))

    return output


def remove_edge_cmap(input_cmap, edge_threshold=5, min_threshold=5):
    """This function removes smeared edges from a coordinate mapping. Depending on the
    interpolation method, coordinate mapping slabs within a larger volume can have
    blurred edges where voxels are interpolated with the neighbouring background. The
    function identifies those affected edge voxels by comparing the difference of each
    voxel to its local neighbourhood. Background is assumed to be filled by zeroes and
    identified edge voxels are set to the background value. A voxels is classified as
    edge outlier if its difference to one local neighbour is larger than edge_threshold
    or if its cmap value is below min_threshold in all dimensions.

    Parameters
    ----------
    input_cmap : str
        Filename of 4d coordinate mapping.
    edge_threshold : float, optional
        Maximum difference to neighbouring voxel (in voxel units). The default
        is 5.
    min_threshold : float, optional
        Minimum cmap value in all dimensions (in voxel units). The default is 5.

    Returns
    -------
    None.

    """
    # load input
    cmap = nb.load(input_cmap)
    cmap_array = cmap.get_fdata()

    # get slab coordinates within volume (assumes background filled with zeroes)
    cmap_slab_coords = np.where(cmap_array != 0)

    # number of voxels
    n_vox = len(cmap_slab_coords[0])
    c_vox = 0

    # initialise for printing out loop status
    c_step = 0
    n_step = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    # identify edges
    for n in range(n_vox):
        # current voxel
        i = cmap_slab_coords[0][n]
        j = cmap_slab_coords[1][n]
        k = cmap_slab_coords[2][n]
        t = cmap_slab_coords[3][n]

        # initialise voxel and its neighbours
        cmap_point = cmap_array[i, j, k, t]
        cmap_neighbour = []

        # consider neighour if not at edge of image volume
        if i > 0:
            cmap_neighbour.append(cmap_array[i - 1, j, k, t])

        if j > 0:
            cmap_neighbour.append(cmap_array[i, j - 1, k, t])

        if k > 0:
            cmap_neighbour.append(cmap_array[i, j, k - 1, t])

        if i < np.shape(cmap_array)[0] - 1:
            cmap_neighbour.append(cmap_array[i + 1, j, k, t])

        if j < np.shape(cmap_array)[1] - 1:
            cmap_neighbour.append(cmap_array[i, j + 1, k, t])

        if k < np.shape(cmap_array)[2] - 1:
            cmap_neighbour.append(cmap_array[i, j, k + 1, t])

        if i > 0 and j > 0:
            cmap_neighbour.append(cmap_array[i - 1, j - 1, k, t])

        if i > 0 and j < np.shape(cmap_array)[1] - 1:
            cmap_neighbour.append(cmap_array[i - 1, j + 1, k, t])

        if i > np.shape(cmap_array)[0] - 1 and j > 0:
            cmap_neighbour.append(cmap_array[i + 1, j - 1, k, t])

        if i > np.shape(cmap_array)[0] - 1 and j < np.shape(cmap_array)[1] - 1:
            cmap_neighbour.append(cmap_array[i + 1, j + 1, k, t])

        if i > 0 and k > 0:
            cmap_neighbour.append(cmap_array[i - 1, j, k - 1, t])

        if i > 0 and k < np.shape(cmap_array)[2] - 1:
            cmap_neighbour.append(cmap_array[i - 1, j, k + 1, t])

        if i > np.shape(cmap_array)[0] - 1 and k > 0:
            cmap_neighbour.append(cmap_array[i + 1, j, k - 1, t])

        if i > np.shape(cmap_array)[0] - 1 and j < np.shape(cmap_array)[2] - 1:
            cmap_neighbour.append(cmap_array[i + 1, j, k + 1, t])

        if j > 0 and k > 0:
            cmap_neighbour.append(cmap_array[i, j - 1, k - 1, t])

        if j > 0 and k < np.shape(cmap_array)[2] - 1:
            cmap_neighbour.append(cmap_array[i, j - 1, k + 1, t])

        if j > np.shape(cmap_array)[1] - 1 and k > 0:
            cmap_neighbour.append(cmap_array[i, j + 1, k - 1, t])

        if j > np.shape(cmap_array)[1] - 1 and j < np.shape(cmap_array)[2] - 1:
            cmap_neighbour.append(cmap_array[i, j + 1, k + 1, t])

        # compare voxel to its neighbours
        cmap_temp = np.abs(np.array(cmap_neighbour) - cmap_point)
        cmap_temp = cmap_temp[~np.isnan(cmap_temp)]
        cmap_temp = np.any(cmap_temp > edge_threshold)

        if cmap_temp:
            cmap_array[i, j, k, t] = np.nan

        # print status
        counter = np.floor(c_vox / n_vox * 100).astype(int)
        if counter == n_step[c_step]:
            print("cmap edge removal: " + str(counter) + " %")
            c_step += 1

        # counter loop cycles
        c_vox += 1

    # remove edges
    cmap_array[np.isnan(cmap_array)] = 0

    # get binary mask from single dimensions
    mask_array = (
        cmap_array[:, :, :, 0] * cmap_array[:, :, :, 1] * cmap_array[:, :, :, 2]
    )
    mask_array[mask_array != 0] = 1

    # get binary mask from min threshold
    min_array1 = cmap_array[:, :, :, 0].copy()
    min_array2 = cmap_array[:, :, :, 1].copy()
    min_array3 = cmap_array[:, :, :, 2].copy()

    min_array1[min_array1 < min_threshold] = 0
    min_array2[min_array2 < min_threshold] = 0
    min_array3[min_array3 < min_threshold] = 0

    min_array1[min_array1 != 0] = 1
    min_array2[min_array2 != 0] = 1
    min_array3[min_array3 != 0] = 1

    min_array = min_array1 + min_array2 + min_array3
    min_array[min_array != 0] = 1

    # mask cmap
    cmap_array[:, :, :, 0] = cmap_array[:, :, :, 0] * mask_array
    cmap_array[:, :, :, 1] = cmap_array[:, :, :, 1] * mask_array
    cmap_array[:, :, :, 2] = cmap_array[:, :, :, 2] * mask_array

    cmap_array[:, :, :, 0] = cmap_array[:, :, :, 0] * min_array
    cmap_array[:, :, :, 1] = cmap_array[:, :, :, 1] * min_array
    cmap_array[:, :, :, 2] = cmap_array[:, :, :, 2] * min_array

    # write output
    path_output, basename_output, ext_output = get_filename(input_cmap)
    output = nb.Nifti1Image(cmap_array, cmap.affine, cmap.header)
    nb.save(output, os.path.join(path_output, basename_output + "_edge" + ext_output))
