# -*- coding: utf-8 -*-

# external inputs
import numpy as np


def get_meshlines(vtx_pial, vtx_white):
    """ Get meshlines
    
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

    Notes
    -------
    created by Daniel Haenelt
    Date created: 24-10-2020
    Last modified: 24-10-2020

    """
    
    # face of line
    fac_new = [[0,1,0]]

    vtx_res = []
    fac_res = []
    for i in range(len(vtx_white)):
    
        # add vertices
        vtx_new = [list(vtx_white[i]),list(vtx_pial[i])]
        
        # update resulting vertex and face list
        vtx_res.extend(vtx_new)
        fac_res.extend(fac_new)
    
        fac_new[0] = [x+2 for x in fac_new[0]]
    
    # vertices and faces as array
    vtx_res = np.array(vtx_res)
    fac_res = np.array(fac_res)
    
    return vtx_res, fac_res
