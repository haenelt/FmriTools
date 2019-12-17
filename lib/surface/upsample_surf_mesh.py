def upsample_surf_mesh(file_in, file_out, niter, method):
    """
    The scripts takes generated FreeSurfer surfaces and upsamples them using Jon Polimeni's function 
    mris_mesh_subdivide. A textfile with additional information (number of vertices and average edge 
    length) is also saved in the same folder.
    Inputs:
        *file_in: filename of input surface.
        *file_out: filename of output surface.
        *niter: number of upsampling iterations.
        *method: upsampling method (linear, loop, butterfly).
    Outputs:
        *orig_params: number of vertices and edge length of original surface.
        *dense_params: number of vertices and edge length of dense surface.

    created by Daniel Haenelt
    Date created: 01-11-2018             
    Last modified: 17-12-2019
    """
    import os
    import numpy as np
    from nibabel.freesurfer.io import read_geometry
    from cortex.polyutils import Surface

    # make output folder
    path_output = os.path.dirname(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # subdivide the surfaces
    os.system("mris_mesh_subdivide" + \
              " --surf " + file_in + \
              " --out " + file_out + \
              " --method " + method + \
              " --iter " + str(niter))

    orig_params = []
    dense_params = []
    pts, polys = read_geometry(file_in)
    pts_dense, polys_dense = read_geometry(file_out)

    # save output parameters (number of vertices and average edge length)
    orig_params.append([np.size(pts[:,0]),
                        Surface(pts,polys).avg_edge_length])

    dense_params.append([np.size(pts_dense[:,0]),
                         Surface(pts_dense,polys_dense).avg_edge_length])
    
    return orig_params, dense_params