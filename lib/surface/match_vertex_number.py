def match_vertex_number(input_white_surf, input_pial_surf, input_white_ind, input_pial_ind):
    """
    This function takes a white and a pial surface as input and matches the vertex numbers of both
    surfaces. Removed vertices (not found in the corresponding other surface) are removed and faces
    are updated. The index files are updated as well.
    Inputs:
        *input_white_surf: filename of white surface.
        *input_pial_surf: filename of pial surface.
        *input_white_ind: corresponding index file of white surface.
        *input_pial_ind: corresponding index file of pial surface.
        
    created by Daniel Haenelt
    Date created: 10-12-2019
    Last modified: 10-12-2019
    """
    import os
    import numpy as np
    from nibabel.freesurfer.io import read_geometry, write_geometry

    # load geometry
    vtx_white, fac_white = read_geometry(input_white_surf)
    vtx_pial, _ = read_geometry(input_pial_surf)

    # load index file
    white_ind = np.loadtxt(input_white_ind).astype(int)
    pial_ind = np.loadtxt(input_pial_ind).astype(int)

    # vertex indices in orig space which are in neither of both deformed surfaces
    ind_remove = list(set(pial_ind) - set(white_ind))
    ind_remove.extend(set(white_ind) - set(pial_ind))
    ind_remove = np.sort(ind_remove) 

    # sort vertices of white surface
    i = 0
    c_white = np.zeros_like(white_ind)
    while i < np.size(white_ind):
        if np.any(ind_remove == white_ind[i]):
            c_white[i] = 1
    
        i += 1    

    # sort vertices of pial surface
    i = 0
    c_pial = np.zeros_like(pial_ind)
    while i < np.size(pial_ind):
        if np.any(ind_remove == pial_ind[i]):
            c_pial[i] = 1

        i += 1  

    # sort faces
    fac_old = fac_white.copy()
    fac_new = fac_white.copy()
    fac_outlier = np.zeros_like(fac_white)
    c_step = 0
    n_step = [10,20,30,40,50,60,70,80,90,100]
    for i in range(len(ind_remove)):
        check_remove = np.where(ind_remove[i] == white_ind)[0]
        if len(check_remove):
            row, col = np.where(fac_old == check_remove)
            fac_outlier[row,col] = 1 # remember which faces to remove
            fac_temp = fac_new.copy() # update face numbering
            fac_temp[fac_old > check_remove] = -1
            fac_temp[fac_temp != -1] = 0
            fac_new += fac_temp
    
        i += 1
        
        # print status
        counter = np.floor(i / len(ind_remove) * 100).astype(int)
        if counter == n_step[c_step]:
            print("sort faces: "+str(counter)+" %")
            c_step += 1

    # remove outlier faces
    fac_outlier = np.sum(fac_outlier,1)
    fac_new = fac_new[fac_outlier == 0]

    # remove outliers in vertices
    vtx_white = vtx_white[c_white == 0]
    vtx_pial = vtx_pial[c_pial == 0]

    # remove outliers in ind
    white_ind = white_ind[c_white == 0]
    pial_ind = pial_ind[c_pial == 0]

    # write output
    path_output = os.path.dirname(input_white_surf)
    name_output = os.path.basename(input_white_surf)
    write_geometry(os.path.join(path_output,name_output+"_match"),vtx_white,fac_new)

    path_output = os.path.dirname(input_pial_surf)
    name_output = os.path.basename(input_pial_surf)
    write_geometry(os.path.join(path_output,name_output+"_match"),vtx_pial,fac_new)
    
    path_output = os.path.dirname(input_white_ind)
    name_output = os.path.splitext(os.path.basename(input_white_ind))[0]
    np.savetxt(os.path.join(path_output,name_output+"_match.txt"), white_ind, fmt='%d')

    path_output = os.path.dirname(input_pial_ind)
    name_output = os.path.splitext(os.path.basename(input_pial_ind))[0]
    np.savetxt(os.path.join(path_output,name_output+"_match.txt"), pial_ind, fmt='%d')