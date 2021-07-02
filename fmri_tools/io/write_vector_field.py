# -*- coding: utf-8 -*-

# external inputs
import numpy as np
from nibabel.freesurfer.io import write_geometry
from gbb.neighbor.nn_2d import nn_2d


def write_vector_field(vtx0, vtx1, fac, adjm, file_out, step_size=100,
                       shape="line"):
    """Write vector field.

    This function generates a surface mesh to visualize a vector field as a 
    triangular mesh.    

    Parameters
    ----------
    vtx0 : ndarray
        Array of vector start points.
    vtx1 : ndarray
        Array of vector end points.
    fac : ndarray
        Corresponding face array.
    adjm : ndarray
        Adjacency matrix.
    file_out : str
        Filename of output surface mesh.
    step_size : int, optional
        Vector subset which will be visualized. The default is 100.
    shape : str, optional
        Line, triangle, prism. The default is "line".

    Returns
    -------
    None.
    
    """

    # array containing a list of considered vectors
    t = np.arange(0, len(vtx0), step_size)

    # initialise faces for specific shape
    if shape == "prism":
        f_new = [[0, 1, 2],
                 [3, 4, 5],
                 [0, 1, 4],
                 [0, 3, 4],
                 [1, 2, 5],
                 [1, 4, 5],
                 [0, 2, 5],
                 [0, 3, 5]]
        f_iter = 6
    elif shape == "triangle":
        f_new = [[0, 1, 2]]
        f_iter = 3
    elif shape == "line":
        f_new = [[0, 1, 0]]
        f_iter = 2

    v_res = []
    f_res = []
    for i in range(len(t)):

        # get index from nearest neighbour of a given vertex
        nn = nn_2d(t[i], adjm, 0)
        nn = nn[:2]

        # get all vertex points for specific shape
        if shape == "prism":
            A = list(vtx0[t[i]])
            B = list(vtx0[nn[0]])
            C = list(vtx0[nn[1]])
            D = list(vtx1[t[i]])
            E = list(vtx1[nn[0]])
            F = list(vtx1[nn[1]])
            v_new = [A, B, C, D, E, F]
        elif shape == "triangle":
            A = list(vtx0[t[i]])
            B = list(vtx0[nn[0]])
            C = list(vtx1[t[i]])
            v_new = [A, B, C]
        elif shape == "line":
            A = list(vtx0[t[i]])
            B = list(vtx1[t[i]])
            v_new = [A, B]

        # update faces
        if i > 0:
            for j in range(len(f_new)):
                f_new[j] = [x + f_iter for x in f_new[j]]

        # update resulting vertex and face list
        v_res.extend(v_new)
        f_res.extend(f_new)

    # vertices and faces as array
    v_res = np.array(v_res)
    f_res = np.array(f_res)

    # write output geometry
    write_geometry(file_out, v_res, f_res)
