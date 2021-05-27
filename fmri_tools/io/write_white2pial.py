# -*- coding: utf-8 -*-

# external inputs
import numpy as np
from nibabel.freesurfer.io import read_geometry, write_geometry
from gbb.neighbor import nn_2d


def write_white2pial(input_white, input_pial, adjm, step_size=100,
                     shape="line"):
    """ Plot white to pial

    This function generates lines between corresponding vertices at the white
    and pial surface to visualize the shift between matched vertices caused by
    realigning surfaces independently. You can either construct prisms,
    triangles or lines.

    Parameters
    ----------
    input_white : str
        Filename of white surface.
    input_pial : str
        Filename of pial surface.
    adjm : obj
        Adjacency matrix.
    step_size : int, optional
        Subset of vertices. The default is 100.
    shape : str, optional
        line, triangle, prism. The default is "line".

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 07-11-2019
    Last modified: 05-10-2020

    """

    # read geometry
    vtx_white, fac_white, header_white = read_geometry(input_white,
                                                       read_metadata=True)
    vtx_pial, fac_pial = read_geometry(input_pial)

    # array containing a list of considered vertices
    t = np.arange(0, len(vtx_white), step_size)

    # initialise faces for specific shape
    if shape == "prism":
        fac_new = [[0, 1, 2],
                   [3, 4, 5],
                   [0, 1, 4],
                   [0, 3, 4],
                   [1, 2, 5],
                   [1, 4, 5],
                   [0, 2, 5],
                   [0, 3, 5]]
        fac_iter = 6
    elif shape == "triangle":
        fac_new = [[0, 1, 2]]
        fac_iter = 3
    elif shape == "line":
        fac_new = [[0, 1, 0]]
        fac_iter = 2

    vtx_res = []
    fac_res = []
    for i in range(len(t)):

        # get index from nearest neighbour of a given vertex
        nn = nn_2d(t[i], adjm, 0)
        nn = nn[:2]

        # get all vertex points for specific shape
        if shape == "prism":
            A = list(vtx_white[t[i]])
            B = list(vtx_white[nn[0]])
            C = list(vtx_white[nn[1]])
            D = list(vtx_pial[t[i]])
            E = list(vtx_pial[nn[0]])
            F = list(vtx_pial[nn[1]])
            vtx_new = [A, B, C, D, E, F]
        elif shape == "triangle":
            A = list(vtx_white[t[i]])
            B = list(vtx_white[nn[0]])
            C = list(vtx_pial[t[i]])
            vtx_new = [A, B, C]
        elif shape == "line":
            A = list(vtx_white[t[i]])
            B = list(vtx_pial[t[i]])
            vtx_new = [A, B]

        # update faces
        if i > 0:
            for j in range(len(fac_new)):
                fac_new[j] = [x + fac_iter for x in fac_new[j]]

        # update resulting vertex and face list
        vtx_res.extend(vtx_new)
        fac_res.extend(fac_new)

    # vertices and faces as array
    vtx_res = np.array(vtx_res)
    fac_res = np.array(fac_res)

    # write output geometry
    write_geometry(input_white + "_plot_white2pial", vtx_res, fac_res,
                   volume_info=header_white)
