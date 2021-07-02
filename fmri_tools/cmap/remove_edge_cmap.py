# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb

# local inputs
from ..io.get_filename import get_filename


def remove_edge_cmap(input_cmap, edge_threshold=5, min_threshold=5):
    """Remove edge cmap.
    
    This function removes smeared edges from a coordinate mapping. Depending on 
    the interpolation method, coordinate mapping slabs within a larger volume 
    can have blurred edges where voxels are interpolated with the neighbouring 
    background. The function identifies those affected edge voxels by comparing 
    the difference of each voxel to its local neighbourhood. Background is 
    assumed to be filled by zeroes and identified edge voxels are set to the 
    background value. A voxels is classified as edge outlier if its difference 
    to one local neighbour is larger than edge_threshold or if its cmap value is 
    below min_threshold in all dimensions.    

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
    mask_array = cmap_array[:, :, :, 0] * cmap_array[:, :, :, 1] * cmap_array[:, :, :, 2]
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
    nb.save(output,
            os.path.join(path_output, basename_output + "_edge" + ext_output))
