# -*- coding: utf-8 -*-

import os

import numpy as np
from nibabel.freesurfer.io import read_geometry, write_geometry


def remove_vertex_outliers(input_surf, input_ind, n=5, overwrite=True):
    """Remove vertex outliers.

    This function removes outlier vertices from a deformed surface mesh. Due
    to interpolation of the deformation fields, edge effects can move some
    vertices to the edge of the image volume. These are removed by comparing
    each vertex to the geometric center of the whole mesh and setting a global
    threshold. The threshold is defined as multiple (n) of the standard
    deviation of the distance distribution. This is a very crude method to
    delete outlier vertices.

    Parameters
    ----------
    input_surf : str
        Path to the input surface mesh.
    input_ind : str
        Path to the corresponding index list (with .txt extension).
    n : TYPE, float
        Threshold parameter. The default is 5.
    overwrite : bool, optional
        Overwrite input surface. The default is True.

    Returns
    -------
    None.

    """

    # load geometry
    vtx, fac = read_geometry(input_surf)

    # load index file
    ind = np.loadtxt(input_ind)

    # get geometry center
    x_mean = np.sum(vtx[:, 0]) / len(vtx[:, 0])
    y_mean = np.sum(vtx[:, 1]) / len(vtx[:, 1])
    z_mean = np.sum(vtx[:, 2]) / len(vtx[:, 2])

    # euclidean distance to geometric center
    vtx_dist = np.zeros_like(vtx)
    vtx_dist[:, 0] = (vtx[:, 0] - x_mean) ** 2
    vtx_dist[:, 1] = (vtx[:, 1] - y_mean) ** 2
    vtx_dist[:, 2] = (vtx[:, 2] - z_mean) ** 2
    vtx_dist = np.sqrt(np.sum(vtx_dist, 1))

    # mean and std distance
    vtx_dist_mean = np.mean(vtx_dist)
    vtx_dist_std = np.std(vtx_dist)

    # distance threshold
    vtx_dist_threshold = vtx_dist_mean + n * vtx_dist_std

    # sort faces
    fac_old = fac.copy()
    fac_outlier = np.zeros_like(fac)
    n_outlier = np.zeros(len(vtx))
    c_step = 0
    n_step = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    for i in range(len(vtx)):
        if vtx_dist[i] > vtx_dist_threshold:
            row, col = np.where(fac_old == i)
            fac_outlier[row, col] = 1  # remember which faces to remove
            n_outlier[i] = 1  # remember which vertices to remove
            fac_temp = fac.copy()  # update face numbering
            fac_temp[fac_old >= i] = -1
            fac_temp[fac_temp != -1] = 0
            fac += fac_temp

        # print status
        counter = np.floor(i / len(vtx) * 100).astype(int)
        if counter == n_step[c_step]:
            print("remove outliers: " + str(counter) + " %")
            c_step += 1

    # remove outlier faces
    fac_outlier = np.sum(fac_outlier, 1)
    fac = fac[fac_outlier == 0]

    # remove outliers in vertex and ind
    vtx = vtx[n_outlier == 0]
    ind = ind[n_outlier == 0]

    # write output
    if overwrite:
        write_geometry(input_surf, vtx, fac)
        np.savetxt(input_ind, ind, fmt="%d")
    else:
        path_output = os.path.dirname(input_surf)
        name_output = os.path.basename(input_surf)
        write_geometry(os.path.join(path_output, name_output + "_out"), vtx, fac)

        path_output = os.path.dirname(input_ind)
        name_output = os.path.splitext(os.path.basename(input_ind))[0]
        np.savetxt(os.path.join(path_output, name_output + "_out.txt"), ind, fmt="%d")
