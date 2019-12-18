def make_mesh(boundary_in, ref_in, file_out, nlayer, niter_smooth=2, niter_upsample=0, 
              niter_inflate=15)
    """
    This function generates a surface mesh from a levelset image. The surface mesh is smoothed and a
    curvature file is generated. Vertices are in the vertex ras coordinate system. Optionally, the
    mesh can be upsampled and an inflated version of the mesh can be written out. The hemisphere
    has to be indicated as prefix in the output file.
    Inputs:
        *boundary_in: filename of 4D levelset image.
        *ref_in: filename of reference volume for getting the coordinate transformation.
        *file_out: filename of output surface.
        *nlayer: layer from the 4D boundary input at which the mesh is generated.
        *niter_smooth: number of smoothing iterations.
        *niter_upsample: number of upsampling iterations (is performed if set > 0).
        *niter_inflate: number of inflating iterations (is performed if set > 0).
    
    created by Daniel Haenelt
    Date created: 18-12-2019
    Last modified: 18-12-2019
    """
    import os
    import nibabel as nb
    from nibabel.affines import apply_affine
    from nighres.surface import levelset_to_mesh
    from nibabel.freesurfer.io import write_geometry
    from cortex.polyutils import Surface
    from lib.surface.vox2ras import vox2ras
    from lib.surface.smooth_surface import smooth_surface
    from lib.surface.upsample_surf_mesh import upsample_surf_mesh
    from lib.surface.get_curvature import get_curvature
    from lib.surface import inflate_surf_mesh
    
    # get levelset boundary from single layer
    boundary = nb.load(boundary_in)
    boundary.header["dim"][0] = 1
    boundary_array = boundary.get_fdata()
    boundary_array = boundary_array[:,:,:,nlayer]
    boundary = nb.Nifti1Image(boundary_array, boundary.affine, boundary.header)
    
    # make mesh
    surf = levelset_to_mesh(boundary, connectivity="18/6", level=0.0, inclusive=True)
    
    # get vertices and faces
    vtx = surf["result"]["points"]
    fac = surf["result"]["faces"]
    
    # get vox2ras transformation
    vox2ras_tkr, _ = vox2ras(ref_in)
    
    # apply vox2ras to vertices
    vtx = apply_affine(vox2ras_tkr, vtx)
    
    # write mesh
    write_geometry(file_out, vtx, fac)
    
    # smooth surface
    smooth_surface(file_out, file_out, niter_smooth)
    
    # upsample mesh (optionally)
    if niter_upsample != 0:
        upsample_surf_mesh(file_out, file_out, niter_upsample, "linear")
        
    # print number of vertices and average edge length
    print("number of vertices: "+str(len(vtx[:,0])))
    print("average edge length: "+str(Surface(vtx,fac).avg_edge_length))
    
    # get curvature (looks for hemisphere prefix)
    get_curvature(file_out, os.path.dirname(file_out))
    
    # inflate surface (optionally)
    if niter_inflate != 0:
        inflate_surf_mesh(file_out, file_out+"_inflated", niter_inflate)