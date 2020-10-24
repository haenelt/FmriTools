# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
from nibabel.freesurfer.io import read_geometry, write_geometry

    
def match_vertex_number(vtx_white, vtx_pial, fac, ind_white, ind_pial):
    """ Match vertex number
    
    This function takes a white and a pial surface as input and matches the 
    vertex numbers of both surfaces. Removed vertices (not found in the 
    corresponding other surface) are removed and faces are updated. The index 
    files are updated as well.    

    Parameters
    ----------
    input_white_surf : str
        Filename of white surface.
    input_pial_surf : str
        Filename of pial surface.
    input_white_ind : str
        Corresponding index file of white surface.
    input_pial_ind : str
        Corresponding index file of pial surface.
    path_output : str
        Path where output is written.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 10-12-2019
    Last modified: 23-10-2020

    """
    
    # vertex indices in orig space which are in neither of both deformed 
    # surfaces
    ind_remove = list(set(ind_pial) - set(ind_white))
    ind_remove.extend(set(ind_white) - set(ind_pial))
    ind_remove = np.sort(ind_remove, rever=True) 

    ind_keep = list(set(fac) - set(ind_remove)).sort()
    
    vtx_white = vtx_white[ind_keep,:]
    vtx_pial = vtx_pial[ind_keep,:]

    # get new faces
    fac_keep = np.zeros(len(fac))
    fac_keep += np.in1d(fac[:,0], ind_keep)
    fac_keep += np.in1d(fac[:,1], ind_keep)
    fac_keep += np.in1d(fac[:,2], ind_keep)
    fac = fac[fac_keep == 3,:]
    
    # reindex faces
    loop_status = 0
    loop_length = len(ind_remove)
    for i in range(loop_length):
        
        # print status
        counter = np.floor(i / loop_length * 100)
        if counter != loop_status:
            print("sort faces: "+str(counter)+" %")
            loop_status = counter
        
        tmp = fac[fac >= ind_remove[i]] - 1
        fac[fac >= ind_remove[i]] = tmp

    return vtx_white, vtx_pial, fac

    # sort vertices of white surface
    #i = 0
    #c_white = np.zeros_like(white_ind)
    #c_pial = np.zeros_like(pial_ind)
    #for i in range(len(ind_remove)):
    #    c_white[white_ind==ind_remove[i]] = 1
    #    c_pial[pial_ind==ind_remove[i]] = 1
    
    #while i < len(white_ind):
    #    if np.any(ind_remove == white_ind[i]):
    #        c_white[i] = 1
    # 
    #    i += 1    

    # sort vertices of pial surface
    #i = 0
    #c_pial = np.zeros_like(pial_ind)
    #while i < len(pial_ind):
    #    if np.any(ind_remove == pial_ind[i]):
    #        c_pial[i] = 1
    #
    #    i += 1  

    # sort faces
    #fac_old = fac_white.copy()
    #fac_new = fac_white.copy()
    #fac_outlier = np.zeros_like(fac_white)
    
    #loop_status = 0
    #loop_length = len(ind_remove)
    #for i in range(loop_length):
    #    check_remove = np.where(ind_remove[i] == white_ind)[0]
    #    if len(check_remove):
    #        row, col = np.where(fac_old == check_remove)
    #        fac_outlier[row,col] = 1 # remember which faces to remove
    #        fac_temp = fac_new.copy() # update face numbering
    #        fac_temp[fac_old > check_remove] = -1
    #        fac_temp[fac_temp != -1] = 0
    #        fac_new += fac_temp
    #
    #    i += 1
    #    
    #    # print status
    #    counter = np.floor(i / loop_length * 100)
    #    if counter != loop_status:
    #        print("sort faces: "+str(counter)+" %")
    #        loop_status = counter

    # remove outlier faces
    #fac_outlier = np.sum(fac_outlier,1)
    #fac_new = fac_new[fac_outlier == 0]

    # remove outliers in vertices
    #vtx_white = vtx_white[c_white == 0]
    #vtx_pial = vtx_pial[c_pial == 0]

    # remove outliers in ind
    #white_ind = white_ind[c_white == 0]
    #pial_ind = pial_ind[c_pial == 0]

    # write output
    #file_out = os.path.join(path_output, os.path.basename(input_white_surf))
    #write_geometry(file_out+"_match", vtx_white, fac_new)

    #file_out = os.path.join(path_output, os.path.basename(input_pial_surf))
    #write_geometry(file_out+"_match", vtx_pial, fac_new)
    
    #file_out = os.path.join(path_output, os.path.basename(input_white_ind))
    #np.savetxt(file_out+"_match.txt", white_ind, fmt='%d')

    #file_out = os.path.join(path_output, os.path.basename(input_pial_ind))
    #np.savetxt(file_out+"_match.txt", pial_ind, fmt='%d')
