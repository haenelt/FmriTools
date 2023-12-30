# -*- coding: utf-8 -*-

import os

import numpy as np
from nibabel.freesurfer.io import read_geometry, write_morph_data
from numpy.linalg import norm
from scipy.stats import sem

from ..io.surf import read_patch


def calculate_distortion(file_patch, file_white, path_output, hemi):
    """Calculate distortion.

    This script computes two metrics (areal distortion and line distortion) to
    estimate the amount of distortion in the flattening process. The form of the
    metric is similar to [1].

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
