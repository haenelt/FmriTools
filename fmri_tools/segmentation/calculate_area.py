# -*- coding: utf-8 -*-

# external inputs
import numpy as np
from numpy.linalg import norm
from nibabel.freesurfer.io import write_morph_data, read_geometry


def calculate_area(filename_surf, filename_area=""):
    """Calculate area.

    The function calculates vertex-wise surface area. The code is taken from the 
    octave code surf2area.m from Anderson Winkler found in his github repository 
    [1]. Consider a triangular face ABC 
    with corner points
    a = [x_A, y_A, z_A]'
    b = [x_B, y_B, z_B]'
    c = [x_C, y_C, z_C]'
    The area for this triangle is given by the normed cross product 
    A = |u x v|/2 with u = a - c and v = b - c. This is a face-wise surface area 
    representation. To convert this to a vertex-wise representation, we assign 
    each vertex one third of the sum of the areas of all faces that meet at that 
    vertex, cf. [2].

    Parameters
    ----------
    filename_surf : str
        Input file geometry on which surface area is calculated.
    filename_area : str, optional
        File name of the surface area file. The default is "".

    Returns
    -------
    dpv : ndarray
        Vertex-wise surface area.

    References
    -------
    .. [1] https://github.com/andersonwinkler/areal
    .. [2] Winkler, A, et al. Measuring and comparing brain cortical surface 
    area and other areal quantities, Neuroimage 61(4), 1428--1443 (2012).
    
    """

    # Read the surface file
    vtx, fac = read_geometry(filename_surf)
    nV = len(vtx)
    nF = len(fac)

    # compute area per face (DPF)
    facvtx = np.concatenate([vtx[fac[:, 0]], vtx[fac[:, 1]], vtx[fac[:, 2]]], axis=1)
    facvtx0 = facvtx[:, 0:6] - np.concatenate([facvtx[:, 6:9], facvtx[:, 6:9]], axis=1)  # place 3rd vtx at origin
    cp = np.cross(facvtx0[:, 0:3], facvtx0[:, 3:6], axisa=1,  axisb=1)  # cross product
    dpf = norm(cp, axis=1) / 2  # half of the norm
    print("Total area (facewise): " + str(np.sum(dpf)))

    # compute area per vertex (DPV)
    dpv = np.zeros(nV)

    # for speed, divide the dpf by 3
    dpf = dpf / 3

    # redistribute
    for f in range(nF):
        dpv[fac[f, :]] = dpv[fac[f, :]] + dpf[f]

    print("Total area (vertexwise): " + str(np.sum(dpv)))

    # save dpv
    if filename_area:
        write_morph_data(filename_area, dpv)

    # return vertex-wise surface area
    return dpv
