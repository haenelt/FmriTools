# -*- coding: utf-8 -*-
"""Layering tools."""

import os

import numpy as np
from cortex.polyutils import Surface
from nibabel.freesurfer.io import read_geometry, write_geometry

from ..surface.mesh import calculate_area

__all__ = ["calc_equivol_surf", "get_meshlines"]


def calc_equivol_surf(file_white, file_pial, n_surfs, factor, niter, hemi, path_output):
    """The script calculates intracortical surfaces based on equi-volumetric layering.
    It is an adaption of Konrad Wagstyl's function in surface_tools. Here, the io_mesh
    is not used anymore and the call to a freesurfer function is omitted. Instead,
    vertex-wise area is calculated in a separate function and we use the nibabel to read
    the surface geometry. First, vertex-wise area is calculated from both input
    geometries. Smoothing to the areas is optional and done if factor is set to a
    non-zero value. Then, based on vertex-wise area, equi-volumetric surfaces are
    computed.

    Parameters
    ----------
    file_white : str
        Input of GM/WM surface.
    file_pial : str
        Input of GM/CSF surface.
    n_surfs : int
        Number of output surfaces (returns input surfaces as 0 and 1).
    factor : float
        Amount of smoothing.
    niter : int
        Number of smoothing iterations.
    hemi : str
        Declare hemisphere for output file.
    path_output : str
        Path where output is saved.
    """
    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # load geometry and area data
    wm_vtx, wm_fac = read_geometry(file_white)
    pial_vtx, pial_fac = read_geometry(file_pial)
    wm_vertexareas = calculate_area(file_white)
    pial_vertexareas = calculate_area(file_pial)

    # smoothing area files (optional)
    if factor != 0:
        wm_vertexareas = Surface(wm_vtx, wm_fac).smooth(
            wm_vertexareas, factor=factor, iterations=niter
        )
        pial_vertexareas = Surface(pial_vtx, pial_fac).smooth(
            pial_vertexareas, factor=factor, iterations=niter
        )

    # number of equally space intracortical surfaces
    vectors = wm_vtx - pial_vtx
    tmp_vtx = pial_vtx.copy()
    tmp_fac = pial_fac.copy()
    mask = vectors.sum(axis=1) != 0  # create mask where vertex coordinates match

    for depth in range(n_surfs):
        print("creating surface " + str(depth + 1))
        betas = _beta(
            float(depth) / (n_surfs - 1), wm_vertexareas[mask], pial_vertexareas[mask]
        )
        betas = np.nan_to_num(betas)
        tmp_vtx[mask] = pial_vtx[mask] + vectors[mask] * np.array([betas]).T
        write_geometry(
            os.path.join(path_output, hemi + "." + "layer" + str(depth)),
            tmp_vtx,
            tmp_fac,
        )


def get_meshlines(vtx_pial, vtx_white):
    """Get meshlines.

    This function returns a vertex and a corresponding face array to visualize
    point-to-point connections between two congruent surface meshs.

    Parameters
    ----------
    vtx_pial : ndarray
        Vertex array of pial surface.
    vtx_white : ndarray
        Vertex array of white surface.

    Returns
    -------
    vtx_res : ndarray
        Vertex array of mesh lines.
    fac_res : ndarray
        Corresponding face array.

    """
    # face of line
    fac_new = [[0, 1, 0]]

    vtx_res = []
    fac_res = []
    for i, _ in enumerate(vtx_white):
        # add vertices
        vtx_new = [list(vtx_white[i]), list(vtx_pial[i])]

        # update resulting vertex and face list
        vtx_res.extend(vtx_new)
        fac_res.extend(fac_new)

        fac_new[0] = [x + 2 for x in fac_new[0]]

    # vertices and faces as array
    vtx_res = np.array(vtx_res)
    fac_res = np.array(fac_res)

    return vtx_res, fac_res


def _beta(alpha, aw, ap):
    """Compute euclidean distance fraction, beta, that will yield the desired volume
    fraction, alpha, given vertex areas in the white matter surface, aw, and on the pial
    surface, ap.

    A surface with `alpha` fraction of the cortical volume below it and `1 - alpha`
    fraction above it can then be constructed from pial, px, and white matter, pw,
    surface coordinates as `beta * px + (1 - beta) * pw`.

    Parameters
    ----------
    alpha : float
        Desired volume fraction.
    aw : ndarray
        Vertex areas in white matter surface.
    ap : ndarray
        Vertex areas in pial surface.

    Returns
    -------
    ndarray
        Euclidean distance fraction.
    """
    if alpha == 0:
        return np.zeros_like(aw)
    if alpha == 1:
        return np.ones_like(aw)

    return 1 - (
        1 / (ap - aw) * (-aw + np.sqrt((1 - alpha) * ap**2 + alpha * aw**2))
    )
