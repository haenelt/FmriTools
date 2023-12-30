# -*- coding: utf-8 -*-

import os

import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
from descartes import PolygonPatch
from numpy.linalg import norm
from scipy.interpolate import griddata
from shapely.geometry import mapping
from skimage.draw import polygon

from ..io.surf import read_patch
from ..segmentation.alpha_shape import alpha_shape


def orthographic_projection(file_patch, img_res, theta, alpha, buffer, path_output):
    """Orthographic projection.

    This script computes a regular grid representation of a flattened patch. It
    is similar to the approach by Kendrick Kay (cvnlookupimages). First, a patch
    is read and the patch coordinates are transformed in relation to its
    barycentre. Optional, a rotation around the barycenter is performed. The
    size of the regular grid is taken from the patch size and is defined with a
    chosen image resolution. The x-axis of the regular grid is flipped to be
    consistent with the RAS coordinate system. Each vertex index is interpolated
    onto the regular grid using nearest neighbours interpolation. A concave hull
    is computed to mask the patch on the regular grid.

    Parameters
    ----------
    file_patch : str
        Filename of flattened patch.
    img_res : float
        Isotropic image resolution in mm.
    theta : float
        Rotation of flat image in deg.
    alpha : float
        Alpha shape value for concave hull computation.
    buffer : float
        Smooth out concave hull.
    path_output : str
        Path where output is saved.

    Returns
    -------
    n_voxel : int
        Number of voxels representing the patch on the regular grid.
    ind_ratio : float
        Ratio of unique indices on the patch.

    """

    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # read patch
    x, y, _, ind = read_patch(file_patch)

    # compute barycentre
    xc = np.sum(x) / np.size(x)
    yc = np.sum(y) / np.size(y)

    # new origin of the patch as vertex with minimum distance to the barycentre
    dist = [norm(np.array((x, y))[:, i] - [xc, yc]) for i in range(len(x))]
    x = x - x[np.argmin(dist)]
    y = y - y[np.argmin(dist)]

    # compute rotation
    theta = np.radians(theta)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c, -s), (s, c)))
    for i in range(len(x)):
        x[i], y[i] = np.dot(R, [x[i], y[i]])

    # target grid to interpolate to
    x_min = np.floor(np.min(x))
    x_max = np.ceil(np.max(x))
    y_min = np.floor(np.min(y))
    y_max = np.ceil(np.max(y))

    xf = np.arange(
        x_max, x_min - img_res, -img_res
    )  # flip x plane to be consistent with RAS
    yf = np.arange(y_min, y_max + img_res, img_res)
    x_plane, y_plane = np.meshgrid(xf, yf)

    # interpolate index values to grid
    x_plane_reshape = x_plane.reshape(
        len(xf) * len(yf),
    )
    y_plane_reshape = y_plane.reshape(
        len(xf) * len(yf),
    )

    # grid interpolation of index data
    coord_orig = np.transpose(np.array((x, y)))
    coord_plane = np.transpose(np.array((x_plane_reshape, y_plane_reshape)))
    method = "nearest"  # can also be linear or cubic
    ind_plane = griddata(coord_orig, ind, coord_plane, method)
    ind_plane = ind_plane.reshape(len(yf), len(xf))

    # get concave hull (alpha shape)
    concave_hull, _ = alpha_shape(coord_orig.tolist(), alpha=alpha)

    # get coordinates of the concave hull
    concave_mapping = mapping(concave_hull.buffer(buffer))
    coord_hull = np.squeeze(np.asarray(concave_mapping["coordinates"]))

    # plot points (swap axes)
    fig1 = plt.figure(figsize=(10, 10))
    ax = fig1.add_subplot(111)
    margin = 0.3
    ax.set_xlim([x_max + margin, x_min - margin])
    ax.set_ylim([y_max + margin, y_min - margin])
    ax.set_xlabel("x in mm")
    ax.set_ylabel("y in mm")
    ax.set_title("Point cloud of flat patch")
    plt.plot(x, y, "o", color="#f16824", markersize=0.5)
    plt.savefig(os.path.join(path_output, os.path.basename(file_patch) + ".points.png"))

    # plot concave hull (swap axes)
    fig2 = plt.figure(
        figsize=(
            10,
            10,
        )
    )
    ax = fig2.add_subplot(111)
    margin = 0.3
    ax.set_xlim([x_max + margin, x_min - margin])
    ax.set_ylim([y_max + margin, y_min - margin])
    ax.set_xlabel("x in mm")
    ax.set_ylabel("y in mm")
    ax.set_title("Concave hull of flat patch")
    polyg = concave_hull.buffer(buffer)
    patch = PolygonPatch(polyg, fc="#999999", ec="#000000", fill=True, zorder=-1)
    ax.add_patch(patch)
    plt.savefig(
        os.path.join(path_output, os.path.basename(file_patch) + ".concave_hull.png")
    )

    # get nearest neighbour coordinates of concave hull on regular grid
    coord_nn = np.zeros_like(coord_hull)
    for i in range(np.size(coord_hull[:, 0])):
        for j in range(np.size(coord_hull[0, :])):
            temp = np.mod(coord_hull[i, j], img_res)
            if temp < img_res / 2:
                coord_nn[i, j] = coord_hull[i, j] - temp
            else:
                coord_nn[i, j] = coord_hull[i, j] + img_res - temp

            if j == 0:
                coord_nn[i, j] = np.argmin(np.abs(xf - coord_nn[i, j]))
            else:
                coord_nn[i, j] = np.argmin(np.abs(yf - coord_nn[i, j]))

    # convert to integer
    coord_nn = coord_nn.astype(int)

    # mask
    mask_plane = np.zeros_like(ind_plane)
    rr, cc = polygon(coord_nn[:, 1], coord_nn[:, 0])
    mask_plane[rr, cc] = 1

    # write niftis
    empty_header = nb.Nifti1Header()
    empty_affine = np.eye(4)

    # new array orientation for niftis
    mask_plane = np.swapaxes(mask_plane, 1, 0)
    ind_plane = np.swapaxes(ind_plane, 1, 0)

    mask = nb.Nifti1Image(mask_plane, empty_affine, empty_header)
    nb.save(mask, os.path.join(path_output, os.path.basename(file_patch) + ".mask.nii"))

    cmap_plane = ind_plane * mask_plane
    cmap = nb.Nifti1Image(cmap_plane, empty_affine, empty_header)
    nb.save(cmap, os.path.join(path_output, os.path.basename(file_patch) + ".cmap.nii"))

    # number of double indices in flat map
    n_voxel = len(cmap_plane[cmap_plane != 0])
    ind_ratio = (len(np.unique(cmap_plane)) - 1) / n_voxel

    return n_voxel, ind_ratio
