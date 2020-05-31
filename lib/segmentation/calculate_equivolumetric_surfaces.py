def calculate_equivolumetric_surfaces(file_white, file_pial, n_surfs, factor, niter, hemi, 
                                      path_output):
    """
    The script calculates intracortical surfaces based on equi-volumetric layering. It is an 
    adaption of Konrad Wagstyl's function in surface_tools. Here, the io_mesh is not used anymore 
    and the call to a freesurfer function is omitted. Instead, vertex-wise area is calculated in a 
    separate function and we use the nibabel to read the surface geometry. First, vertex-wise area 
    is calculated from both input geometries. Smoothing to the areas is optional and done if factor 
    is set to a non-zero value. Then, based on vertex-wise area, equi-volumetric surfaces are 
    computed.
    Inputs:
        *file_white: input of GM/WM surface.
        *file_pial: input of GM/CSF surface.
        *n_surfs: number of output surfaces (returns input surfaces as 0 and 1).
        *factor: amount of smoothing.
        *niter: number of smoothing iterations.
        *hemi: declare hemisphere for output file.
        *path_output: path where output is saved.
    
    created by Daniel Haenelt
    Date created: 01-11-2018             
    Last modified: 17-12-2018
    """
    import os
    import numpy as np
    from nibabel.freesurfer.io import read_geometry, write_geometry 
    from cortex.polyutils import Surface
    from lib.segmentation.calculate_area import calculate_area

    def beta(alpha, aw, ap):
        """Compute euclidean distance fraction, beta, that will yield the desired
        volume fraction, alpha, given vertex areas in the white matter surface, 
        aw, and on the pial surface, ap.
    
        A surface with `alpha` fraction of the cortical volume below it and 
        `1 - alpha` fraction above it can then be constructed from pial, px, and 
        white matter, pw, surface coordinates as `beta * px + (1 - beta) * pw`.
        """
        if alpha == 0:
            return np.zeros_like(aw)
        elif alpha == 1:
            return np.ones_like(aw)
        else:
            return 1-(1 / (ap - aw) * (-aw + np.sqrt((1-alpha)*ap**2 + alpha*aw**2)))

    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # load geometry and area data
    wm_vtx, wm_fac = read_geometry(file_white)
    pial_vtx, pial_fac = read_geometry(file_pial)
    wm_vertexareas = calculate_area(file_white)
    pial_vertexareas = calculate_area(file_pial)

    # smoothing area files (optional)
    if factor != 0:
        wm_vertexareas = Surface(wm_vtx,wm_fac).smooth(wm_vertexareas, factor=factor, iterations=niter)
        pial_vertexareas = Surface(pial_vtx,pial_fac).smooth(pial_vertexareas, factor=factor, iterations=niter)

    # number of equally space intracortical surfaces
    vectors = wm_vtx - pial_vtx
    tmp_vtx = pial_vtx.copy()
    tmp_fac = pial_fac.copy()
    mask = vectors.sum(axis=1) != 0 # create mask where vertex coordinates match
    
    for depth in range(n_surfs):
        print("creating surface " + str(depth +1))
        betas = beta(float(depth)/(n_surfs-1), wm_vertexareas[mask], pial_vertexareas[mask])
        betas = np.nan_to_num(betas)
        tmp_vtx[mask] = pial_vtx[mask] + vectors[mask]* np.array([betas]).T
        write_geometry(os.path.join(path_output,hemi+"."+"layer"+str(depth)), tmp_vtx,tmp_fac)
