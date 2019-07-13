def upsample_surf_mesh(path, sub, hemi, niter, method, path_output):
    """
    The scripts takes generated FreeSurfer surfaces and upsamples them using Jon Polimeni's function 
    mris_mesh_subdivide. A textfile with additional information (number of vertices and average edge 
    length) is also saved in the same folder.
    Inputs:
        *path: path to the freesurfer segmentation folder.
        *sub: name of the freesurfer segmentation folder.
        *hemi: hemisphere.
        *niter: number of upsampling iterations.
        *method: upsampling method (linear, loop, butterfly).
        *path_output: path where output is saved.
    Outputs:
        *orig_params: number of vertices and edge length of original surface.
        *dense_params: number of vertices and edge length of dense surface.

    created by Daniel Haenelt
    Date created: 01-11-2018             
    Last modified: 13-07-2019
    """
    import os
    import numpy as np
    from nibabel.freesurfer.io import read_geometry
    from cortex.polyutils import Surface
  
    # parameters
    surfs = ["sphere", "white", "pial", "inflated"] # list of surfaces to subdivide

    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # get filenames for both hemispheres
    filename = [hemi + "." + surfs for surfs in surfs]

    # subdivide the surfaces
    for i in range(len(filename)):
        os.system("mris_mesh_subdivide" + \
                  " --surf " + os.path.join(path,sub,"surf",filename[i]) + \
                  " --out " + os.path.join(path_output,filename[i]) + \
                  " --method " + method + \
                  " --iter " + str(niter))

    # write information of the upsampled mesh to textfile
    type = "white" # choose white surface

    orig_params = []
    dense_params = []
    pts, polys = read_geometry(os.path.join(path,sub,"surf",hemi+"."+type))
    pts_dense, polys_dense = read_geometry(os.path.join(path_output,hemi+"."+type))

    # save output parameters (number of vertices and average edge length)
    orig_params.append([np.size(pts[:,0]),
                        Surface(pts,polys).avg_edge_length])

    dense_params.append([np.size(pts_dense[:,0]),
                         Surface(pts_dense,polys_dense).avg_edge_length])
    
    return orig_params, dense_params