def mesh_sampling_other(surf_in, file_in, source_in, target_in, target2source_in, source2target_in,
                        boundaries_in, path_output, layer, smooth_iter=2, r=[0.4,0.4,0.4],
                        interpolation="Cu", average_layer=False, write_profile=False, 
                        write_upsampled=True, cleanup=True):
    """
    This function samples data from an image volume to a surface mesh which is located in a 
    different space. Boundaries and surface mesh are first transformed to the space of the image 
    volume using coordinate mappings before data sampling. If average_layer is true, the parameter 
    layer should contain only two integers which denote the start and ending layer. The basename
    of the surface file should have no file extension and the hemisphere should be stated as prefix.
    Inputs:
        *surf_in: filename of input surface mesh.
        *file_in: filename of input volume from which data is sampled.
        *source_in: reference source file.
        *target_in: reference target file.
        *target2source_in: target to source coordinate mapping.
        *source2target_in: source to target coordinate mapping.
        *boundaries_in: filename of 4D levelset image.
        *path_output: path where output is written.
        *layer: which layers to sample (array of integers).
        *smooth_iter: number of smoothing iterations after mesh deformation.
        *r: destination voxel size after upsampling (performed if not None).
        *interpolation: interpolation method for upsampling of file from whic data is sampled.
        *average_layer: average across cortex.
        *write_profile: write sampled profile.
        *write_upsampled: write upsampled file.
        *cleanup: remove intermediate files.
    
    created by Daniel Haenelt
    Date created: 13-01-2020
    Last modified: 14-01-2020
    """
    import os
    import shutil as sh
    import nibabel as nb
    from nighres.registration import apply_coordinate_mappings
    from lib.io.get_filename import get_filename
    from lib.utils.upsample_volume import upsample_volume
    from lib.surface.deform_surface import deform_surface
    from lib.surface.mesh_sampling import mesh_sampling
       
    """
    set folder structure
    """
    path_temp = os.path.join(path_output,"temp")
    path_cmap = os.path.join(path_temp,"cmap")
    path_data = os.path.join(path_temp,"data")
    path_surf = os.path.join(path_temp,"surf")
     
    if not os.path.exists(path_output):
        os.makedirs(path_output)
    
    if not os.path.exists(path_temp):
        os.makedirs(path_temp)
    
    if not os.path.exists(path_cmap):
        os.makedirs(path_cmap)
        
    if not os.path.exists(path_data):
        os.makedirs(path_data)

    if not os.path.exists(path_surf):
        os.makedirs(path_surf)

    # get filenames
    _, hemi, name_surf = get_filename(surf_in)
    _, name_file, ext_file = get_filename(file_in)
    _, name_t2s, ext_t2s = get_filename(target2source_in)
    _, name_s2t, ext_s2t = get_filename(source2target_in)
    _, name_target, ext_target = get_filename(target_in)
    _, name_source, ext_source = get_filename(source_in)
    
    # copy input files
    sh.copyfile(target2source_in, os.path.join(path_cmap, "t2s"+ext_t2s))
    sh.copyfile(source2target_in, os.path.join(path_cmap, "s2t"+ext_s2t))
    sh.copyfile(file_in, os.path.join(path_data, name_file+ext_file))
    sh.copyfile(target_in, os.path.join(path_data, "target"+ext_target))
    sh.copyfile(source_in, os.path.join(path_data, "source"+ext_source))
    
    # set filenames
    data = os.path.join(path_data, name_file+ext_file)
    data_upsampled = os.path.join(path_data, name_file+"_upsampled"+ext_file)
    source = os.path.join(path_data, "source"+ext_source)
    source_upsampled = os.path.join(path_data, "source_upsampled"+ext_source)
    target = os.path.join(path_data, "target"+ext_target)
    target_upsampled = os.path.join(path_data, "target_upsampled"+ext_target)
    t2s = os.path.join(path_cmap, "t2s"+ext_t2s)
    t2s_upsampled = os.path.join(path_cmap, "t2s_upsampled"+ext_t2s)
    t2s_rescaled = os.path.join(path_cmap, "t2s_upsampled_rescaled"+ext_t2s)
    s2t = os.path.join(path_cmap, "s2t"+ext_t2s)
    s2t_upsampled = os.path.join(path_cmap,"s2t_upsampled"+ext_s2t)
    s2t_rescaled = os.path.join(path_cmap, "s2t_upsampled_rescaled"+ext_s2t)
    
    """
    upsample data
    """
    if r:
        upsample_volume(data, data_upsampled, dxyz = r, rmode = interpolation)
        upsample_volume(t2s, t2s_upsampled, dxyz = r, rmode = "Linear")
        upsample_volume(s2t, s2t_upsampled, dxyz = r, rmode = "Linear")
        upsample_volume(target, target_upsampled, dxyz = r, rmode = "Linear")
        upsample_volume(source, source_upsampled, dxyz = r, rmode = "Linear")
    
    """
    rescale cmap
    """
    dim_target = nb.load(target).header["dim"][1:4] - 1
    dim_source = nb.load(source).header["dim"][1:4] - 1
    dim_target_upsampled = nb.load(target_upsampled).header["dim"][1:4] - 1
    dim_source_upsampled = nb.load(source_upsampled).header["dim"][1:4] - 1
    
    cmap_t2s = nb.load(t2s_upsampled)
    cmap_s2t = nb.load(s2t_upsampled)
    cmap_t2s_array = cmap_t2s.get_fdata()
    cmap_s2t_array = cmap_s2t.get_fdata()
    
    for i in range(3):
        cmap_t2s_array[:,:,:,i] = cmap_t2s_array[:,:,:,i] / dim_target[i] * dim_target_upsampled[i]
        cmap_s2t_array[:,:,:,i] = cmap_s2t_array[:,:,:,i] / dim_source[i] * dim_source_upsampled[i]
    
    cmap_t2s_rescaled = nb.Nifti1Image(cmap_t2s_array, cmap_t2s.affine, cmap_t2s.header)
    cmap_s2t_rescaled = nb.Nifti1Image(cmap_s2t_array, cmap_s2t.affine, cmap_s2t.header)
    nb.save(cmap_t2s_rescaled, t2s_rescaled)
    nb.save(cmap_s2t_rescaled, s2t_rescaled)
    
    """
    deform boundaries
    """
    apply_coordinate_mappings(image = boundaries_in, # input 
                              mapping1 = t2s_rescaled, # cmap
                              interpolation = "linear", # nearest or linear
                              padding = "zero", # closest, zero or max
                              save_data = True, # save output data to file (boolean)
                              overwrite = True, # overwrite existing results (boolean)
                              output_dir = path_data, # output directory
                              file_name = "boundaries" # base name with file extension for output
                              )
    
    """
    deform mesh
    """
    deform_surface(input_surf = surf_in, 
                   input_orig = target_upsampled, 
                   input_deform = s2t_rescaled,
                   input_target = source_upsampled, 
                   hemi = hemi,
                   path_output = path_surf,
                   smooth_iter = smooth_iter,
                   sort_faces = False, 
                   cleanup = True)
    
    """
    sample data
    """    
    if smooth_iter:
        surf_in = os.path.join(path_surf, hemi+name_surf+"_def_smooth")
    else:
        surf_in = os.path.join(path_surf, hemi+name_surf+"_def")
    
    mesh_sampling(surf_in = surf_in, 
                  file_in = data_upsampled, 
                  boundaries_in = os.path.join(path_data, "boundaries_def-img.nii.gz"), 
                  path_output = path_output, 
                  layer = layer,
                  r = None, 
                  average_layer = average_layer, 
                  write_profile = write_profile,
                  write_upsampled = write_upsampled)
    
    # delete intermediate files
    if cleanup:
        sh.rmtree(path_temp, ignore_errors=True)
