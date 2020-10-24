# -*- coding: utf-8 -*-

# external inputs
import numpy as np

    
def match_vertex_number(vtx_white, vtx_pial, fac, ind_white, ind_pial):
    """ Match vertex number
    
    This function takes arrays of white and pial surface vertices as input and 
    matches the vertex numbers of both surfaces. The match is based on index
    lists which map the vertex indices to a common reference space. Removed 
    vertices (not found in the corresponding other surface) are removed and the 
    face array is updated. Index files are updated as well.

    Parameters
    ----------
    vtx_white : ndarray
        Vertex array of white surface.
    vtx_pial : ndarray
        Vertex array of pial surface.
    fac : ndarray
        Corresponding face array.
    ind_white : list
        Index list of white surface.
    ind_pial : list
        Index list of pial surface.

    Returns
    -------
    vtx_white : ndarray
        Updated vertex array of white surface.
    vtx_pial : ndarray
        Updated vertex array of pial surface.
    fac_new : ndarray
        Updated face array.
    ind_white : list
        updated index list.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 10-12-2019
    Last modified: 23-10-2020

    """
      
    # vertex indices in common reference space which are in neither of both 
    # deformed surfaces
    ind_remove = list(set(ind_pial) - set(ind_white))
    ind_remove.extend(set(ind_white) - set(ind_pial))
    ind_remove = sorted(ind_remove, reverse=True)

    # sort vertices of white surface
    print("sort white surface vertices")
    
    i = 0
    c_white = np.zeros_like(ind_white)
    while i < len(ind_white):
        if np.any(ind_remove == ind_white[i]):
            c_white[i] = 1
     
        i += 1    

    # sort vertices of pial surface
    print("sort pial surface vertices")
    
    i = 0
    c_pial = np.zeros_like(ind_pial)
    while i < len(ind_pial):
        if np.any(ind_remove == ind_pial[i]):
            c_pial[i] = 1
    
        i += 1  

    # sort faces
    fac_old = fac.copy()
    fac_new = fac.copy()
    fac_outlier = np.zeros_like(fac)
    
    loop_status = 0
    loop_length = len(ind_remove)
    for i in range(loop_length):
        check_remove = np.where(ind_remove[i] == ind_white)[0]
        if len(check_remove):
            row, col = np.where(fac_old == check_remove)
            fac_outlier[row,col] = 1 # remember which faces to remove
            fac_temp = fac_new.copy() # update face numbering
            fac_temp[fac_old > check_remove] = -1
            fac_temp[fac_temp != -1] = 0
            fac_new += fac_temp
    
        i += 1
        
        # print status
        counter = np.floor(i / loop_length * 100)
        if counter != loop_status:
            print("sort faces: "+str(counter)+" %")
            loop_status = counter

    # remove outliers in faces
    fac_outlier = np.sum(fac_outlier,1)
    fac_new = fac_new[fac_outlier == 0]

    # remove outliers in vertices
    vtx_white = vtx_white[c_white == 0]
    vtx_pial = vtx_pial[c_pial == 0]

    # remove outliers in ind
    ind_white = ind_white[c_white == 0]
    
    return vtx_white, vtx_pial, fac_new, ind_white
