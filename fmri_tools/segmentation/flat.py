# -*- coding: utf-8 -*-
"""Flattening utilities."""

import datetime
import math
import os
import shutil as sh

import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
import shapely.geometry as geometry
from descartes import PolygonPatch
from nibabel.freesurfer.io import read_geometry, write_morph_data
from numpy.linalg import norm
from scipy.interpolate import griddata
from scipy.spatial.qhull import Delaunay
from scipy.stats import sem
from shapely.geometry import mapping
from shapely.ops import cascaded_union, polygonize
from skimage.draw import polygon

from .. import execute_command
from ..io.surf import label_as_patch, read_patch

__all__ = [
    "orthographic_projection",
    "alpha_shape",
    "label_flattening",
    "surface_flattening",
    "calculate_distortion",
]


def orthographic_projection(file_patch, img_res, theta, alpha, buffer, path_output):
    """This script computes a regular grid representation of a flattened patch. It is
    similar to the approach by Kendrick Kay (cvnlookupimages). First, a patch is read
    and the patch coordinates are transformed in relation to its barycentre. Optional, a
    rotation around the barycenter is performed. The size of the regular grid is taken
    from the patch size and is defined with a chosen image resolution. The x-axis of the
    regular grid is flipped to be consistent with the RAS coordinate system. Each vertex
    index is interpolated onto the regular grid using nearest neighbours interpolation.
    A concave hull is computed to mask the patch on the regular grid.

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
    for i, _ in enumerate(x):
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


def alpha_shape(points, alpha):
    """This function computes the alpha shape (concave hull) of a set of points. It is
    taken and slightly changed from [1]. Delaunay triangles are computed which establish
    a connection between each point and nearby points and then we remove some of the
    triangles that are too far from their neighbors. This removal part is the key. By
    identifying candidates for removal we are saying that these points are too far from
    their connected points so don't use that connection as part of the boundary.

    Parameters
    ----------
    points : list
        Iterable container of points.
    alpha : float
        Alpha value to influence the gooeyness of the border. Smaller numbers
        don't fall inward as much as larger numbers. Too large, and you lose
        everything!

    Returns
    -------
    cascaded_union(triangles) : poly
        Concave hull coordinates.
    edge_points : list
        Start and end points of all edges.

    References
    -------
    .. [1] http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/
    (accessed 22-10-2018)

    """
    # when you have one triangle, there is no sense in computing an alpha shape
    if len(points) < 4:
        return geometry.MultiPoint(list(points)).convex_hull

    # add a line between point i and j if not alread added
    def add_edge(edges, edge_points, coords, i, j):
        if (i, j) in edges or (j, i) in edges:
            return
        edges.add((i, j))
        edge_points.append(coords[[i, j]])

    # loop over triangles (ia, ib, ic are indices of triangle corner points)
    coords = np.array(points)
    tri = Delaunay(coords)
    edges = set()
    edge_points = []
    for ia, ib, ic in tri.vertices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]

        # lengths of sides of triangle
        a = math.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = math.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = math.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)

        # semiperimeter of triangle
        s = (a + b + c) / 2.0

        # area of triangle by Heron's formula
        area = math.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)

        # radius filter.
        if circum_r < 1.0 / alpha:
            add_edge(edges, edge_points, coords, ia, ib)
            add_edge(edges, edge_points, coords, ib, ic)
            add_edge(edges, edge_points, coords, ic, ia)

    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))

    return cascaded_union(triangles), edge_points


def label_flattening(file_ref, file_label, path_output, cleanup=True):
    """Uses the FreeSurfer mris_flatten function to flatten a dense patch of manually
    defined cortex. A label file is used to define the flattened region. Instead of
    smoothwm, we use white for surface flattening.

    Parameters
    ----------
    file_ref : str
        Reference surface file for flattening.
    file_label : str
        Label to be flattened saved as <hemi>.<namePATCH>.label.
    path_output : str
        Path where output is saved.
    cleanup : bool, optional
        Delete intermediate files. The default is True.
    """
    # create temporary folder
    tmp1 = np.random.randint(0, 10, 5)
    tmp1 = "".join(str(i) for i in tmp1)
    tmp2 = datetime.datetime.now().strftime("%S%f")
    tmp_string = tmp1 + tmp2
    path_temp = os.path.join(path_output, "tmp_" + tmp_string)

    # make temporary folder
    if not os.path.exists(path_temp):
        os.makedirs(path_temp)
    else:
        raise FileExistsError("Temporary folder already exists!")

    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # change to temporary folder
    cwd = os.getcwd()
    os.chdir(path_temp)

    # get hemisphere from file name
    hemi = os.path.basename(file_label)[:2]
    name_patch = "." + os.path.basename(file_label).split(".")[1]

    # convert label to patch file
    file_patch = os.path.join(path_temp, os.path.join(hemi + name_patch + ".patch.3d"))
    label_as_patch(file_ref, file_label, file_patch)

    # copy reference file and path into temporary folder
    sh.copy2(file_ref, os.path.join(path_temp, hemi + ".smoothwm"))

    # surface flattening
    w = 0  # write out the surface every number of iterations.
    s = 20  # size of neighbourhood to be used in the optimization
    n = 7  # number of vertices at each distance to be used in the optimization
    command = "mris_flatten"
    command += f" -w {w}"
    command += f" -distances {s} {n}"
    command += f" {file_patch}"
    command += f" {hemi}{name_patch}.patch.flat"

    # run
    execute_command(command)

    # copy output
    sh.copy2(
        os.path.join(path_temp, hemi + name_patch + ".patch.flat"),
        os.path.join(path_output, hemi + name_patch + ".patch.flat"),
    )
    sh.copy2(
        os.path.join(path_temp, hemi + name_patch + ".patch.flat.out"),
        os.path.join(path_output, hemi + name_patch + ".patch.flat.out"),
    )

    # change to old folder
    os.chdir(cwd)

    # delete temporary files
    if cleanup:
        sh.rmtree(path_temp, ignore_errors=True)


def surface_flattening(file_ref, file_patch, path_output, cleanup=True):
    """Uses the FreeSurfer mris_flatten function to flatten a dense patch of manually
    defined cortex. The manually cutted patch should have the following file name
    <hemi>.<namePATCH>.patch.3d. Instead of smoothwm, we use white for surface
    flattening.

    Parameters
    ----------
    file_ref : str
        Reference surface file for flattening.
    file_patch : str
        Patch to be flattened saved as <hemi>.<namePATCH>.patch.3d.
    path_output : str
        Path where output is saved.
    cleanup : bool, optional
        Delete intermediate files. The default is True.
    """
    # create temporary folder
    tmp1 = np.random.randint(0, 10, 5)
    tmp1 = "".join(str(i) for i in tmp1)
    tmp2 = datetime.datetime.now().strftime("%S%f")
    tmp_string = tmp1 + tmp2
    path_temp = os.path.join(path_output, "tmp_" + tmp_string)

    # make temporary folder
    if not os.path.exists(path_temp):
        os.makedirs(path_temp)
    else:
        raise FileExistsError("Temporary folder already exists!")

    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # change to temporary folder
    cwd = os.getcwd()
    os.chdir(path_temp)

    # divide patch basename
    hemi = os.path.splitext(
        os.path.splitext(os.path.splitext(os.path.basename(file_patch))[0])[0]
    )[0]
    name_patch = os.path.splitext(
        os.path.splitext(os.path.splitext(os.path.basename(file_patch))[0])[0]
    )[1]

    # copy reference file and path into temporary folder
    sh.copy2(file_ref, os.path.join(path_temp, hemi + ".smoothwm"))
    sh.copy2(file_patch, os.path.basename(file_patch))

    # surface flattening
    w = 0  # write out the surface every number of iterations.
    s = 20  # size of neighbourhood to be used in the optimization
    n = 7  # number of vertices at each distance to be used in the optimization
    command = "mris_flatten"
    command += f" -w {w}"
    command += f" -distances {s} {n}"
    command += f" {hemi}{name_patch}.patch.3d"
    command += f" {hemi}{name_patch}.patch.flat"

    # run
    execute_command(command)

    # copy output
    sh.copy2(
        os.path.join(path_temp, hemi + name_patch + ".patch.flat"),
        os.path.join(path_output, hemi + name_patch + ".patch.flat"),
    )
    sh.copy2(
        os.path.join(path_temp, hemi + name_patch + ".patch.flat.out"),
        os.path.join(path_output, hemi + name_patch + ".patch.flat.out"),
    )

    # change to old folder
    os.chdir(cwd)

    # delete temporary files
    if cleanup:
        sh.rmtree(path_temp, ignore_errors=True)


def calculate_distortion(file_patch, file_white, path_output, hemi):
    """This script computes two metrics (areal distortion and line distortion) to
    estimate the amount of distortion in the flattening process. The form of the metric
    is similar to [1].

    For vertex-wise areal distortion (VAD), first all faces which lie within the
    patch are searched. Then, triangle areas before and after flattening are
    computed for each face. The amount of distortion is defined as ratio between
    both areas. A vertex-wise representation is computed by considereing each
    vertex point within the patch and taking the sum of areal distortions of all
    neighbouring faces.

    For vertex-wise line distortion (VLD), for each node in the patch, we take
    the vertex points before and after flattening and look for all neighbouring
    nodes, i.e., we are looking for all faces which contain the node. We take
    then the sum of all euclidean distances to all neighbouring nodes and before
    and after flattening separately and compute the ratio between summed
    distances.

    Note that only points within the patch are taken into account which have
    full faces within the point cloud of the patch corresponding to the original
    white surface. I.e., a few border points are excluded from the analysis and
    set to zero in the morphological output file. The number of excluded points
    is returned for each metric.

    Parameters
    ----------
    file_patch : str
        Filename of flattened patch.
    file_white : str
        Filename of white surface.
    path_output : str
        Path where output is written.
    hemi : str
        Hemisphere.

    Returns
    -------
    VAD_params : list
        Descriptive parameters of areal distortion.
    VLD_params : list
        Descriptive parameters of line distortion.

    References
    -------
    .. [1] http://brainvis.wustl.edu/wiki/index.php/Caret:Operations/Morphing

    """
    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # load data
    vtx_white, fac_white = read_geometry(file_white)
    _, _, _, ind_patch = read_patch(file_patch)
    vtx_patch = np.zeros((len(ind_patch), 3))
    vtx_patch[:, 0], vtx_patch[:, 1], vtx_patch[:, 2], _ = read_patch(file_patch)

    # look for faces which exist in the patch
    fac_patch = []
    for i in range(len(fac_white)):
        if (
            np.any(fac_white[i, 0] == ind_patch)
            and np.any(fac_white[i, 1] == ind_patch)
            and np.any(fac_white[i, 2] == ind_patch)
        ):
            fac_patch.append(fac_white[i, :])
    fac_patch = np.array(fac_patch)

    # Areal distortion

    # calculate face-wise areal distortion (before flattening)
    facvtx_white = np.concatenate(
        [
            vtx_white[fac_patch[:, 0]],
            vtx_white[fac_patch[:, 1]],
            vtx_white[fac_patch[:, 2]],
        ],
        axis=1,
    )
    facvtx0_white = facvtx_white[:, 0:6] - np.concatenate(
        [facvtx_white[:, 6:9], facvtx_white[:, 6:9]], axis=1
    )  # place 3rd vtx at origin
    cp = np.cross(
        facvtx0_white[:, 0:3], facvtx0_white[:, 3:6], axisa=1, axisb=1
    )  # cross product
    A_white = norm(cp, axis=1) / 2  # half of the norm

    # calculate face-wise areal distortion (after flattening)
    vtx_patch_all = np.zeros_like(vtx_white).astype(float)
    for i in range(len(ind_patch)):
        vtx_patch_all[ind_patch[i], :] = vtx_patch[i, :]

    facvtx_patch = np.concatenate(
        [
            vtx_patch_all[fac_patch[:, 0]],
            vtx_patch_all[fac_patch[:, 1]],
            vtx_patch_all[fac_patch[:, 2]],
        ],
        axis=1,
    )
    facvtx0_patch = facvtx_patch[:, 0:6] - np.concatenate(
        [facvtx_patch[:, 6:9], facvtx_patch[:, 6:9]], axis=1
    )  # place 3rd vtx at origin
    cp = np.cross(
        facvtx0_patch[:, 0:3], facvtx0_patch[:, 3:6], axisa=1, axisb=1
    )  # cross product
    A_patch = norm(cp, axis=1) / 2  # half of the norm

    # calculate face-wise distortion
    A_dist = A_patch / A_white

    # convert to vertex-wise representation
    VAD = np.zeros(len(vtx_white)).astype(float)
    VAD_miss = 0
    for i in range(len(ind_patch)):
        temp = np.where(fac_patch == ind_patch[i])
        if np.any(temp[0]):
            VAD[ind_patch[i]] = np.sum(A_dist[temp[0]]) / len(temp[0])
        else:
            VAD_miss += 1

    VAD_params = [
        np.mean(VAD[ind_patch]),
        np.std(VAD[ind_patch]),
        sem(VAD[ind_patch]),
        np.min(VAD[ind_patch]),
        np.max(VAD[ind_patch]),
        VAD_miss,
    ]

    # save morphological data
    write_morph_data(
        os.path.join(path_output, os.path.basename(file_patch) + ".areal_distortion"),
        VAD,
    )

    # Linear distortion

    VLD = np.zeros(len(vtx_white)).astype(float)
    VLD_miss = 0
    for i in range(len(ind_patch)):
        node_patch = vtx_patch_all[ind_patch[i]]
        node_white = vtx_white[ind_patch[i]]

        ind_temp = fac_patch[np.where(fac_patch == ind_patch[i])[0], :]
        ind_temp = np.reshape(ind_temp, len(ind_temp) * 3)
        ind_temp = np.unique(ind_temp)
        ind_temp = np.delete(ind_temp, np.where(ind_temp == ind_patch[i]))

        if np.any(ind_temp):
            VLD_patch = 0
            VLD_white = 0
            for j in range(len(ind_temp)):
                VLD_patch += norm(node_patch - vtx_patch_all[ind_temp[j]])
                VLD_white += norm(node_white - vtx_white[ind_temp[j]])

            VLD[ind_patch[i]] = VLD_patch / VLD_white
        else:
            VLD_miss += 1

    VLD_params = [
        np.mean(VLD[ind_patch]),
        np.std(VLD[ind_patch]),
        sem(VLD[ind_patch]),
        np.min(VLD[ind_patch]),
        np.max(VLD[ind_patch]),
        VLD_miss,
    ]

    # save morphological data
    write_morph_data(
        os.path.join(path_output, os.path.basename(file_patch) + ".line_distortion"),
        VLD,
    )

    return VAD_params, VLD_params
