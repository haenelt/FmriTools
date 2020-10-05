def get_b0_orientation(surf_in, vol_in, write_output=False, path_output="", name_output=""):
    """
    This function computes the angle between surface normals and B0-direction per vertex.
    Inputs:
        *surf_in: input of surface mesh.
        *vol_in: input of corresponding nifti volume.
        *write output: write out to disk (boolean).
        *path_output: path where to save output.
        *name_output: basename of output file.
    Outputs:
        *theta: angle in radians.
        
    created by Daniel Haenelt
    Date created: 31-07-2020 
    Last modified: 05-10-2020
    """
    import os
    import numpy as np
    import nibabel as nb
    from nibabel.affines import apply_affine
    from nibabel.freesurfer.io import read_geometry, write_morph_data
    from lib.io.get_filename import get_filename
    from lib.surface.vox2ras import vox2ras
    from gbb.normal import get_normal
    
    # make subfolders
    if write_output and not os.path.exists(path_output):
        os.makedirs(path_output)
    
    # get hemi from surface filename
    _, hemi, _ = get_filename(surf_in)
    
    # load surface
    vtx, fac = read_geometry(surf_in)
    
    # get transformation matrix
    _, r2v = vox2ras(vol_in)      # ras-tkr -> voxel
    v2s = nb.load(vol_in).affine  # voxel -> scanner-ras
    M = v2s.dot(r2v)
    
    # apply affine transformation
    vtx = apply_affine(M, vtx)
    
    # get surface normals
    n = get_normal(vtx, fac)
    
    # get angle between b0 and surface normals in radians    
    theta = np.arccos(np.dot(n, [0,0,1]))
    
    # write output
    if write_output:
        write_morph_data(os.path.join(path_output, hemi+"."+name_output), theta)

    return theta