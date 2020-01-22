def deform_surface(input_surf, input_orig, input_deform, input_target, hemi, path_output, 
                   smooth_iter=0, sort_faces=False, flip_faces=False, cleanup=True):
    """
    This function deforms a surface mesh in freesurfer convention using a coordinate map containing
    voxel coordinates. The computation takes quite a while because in the case of removed vertices,
    the remaining faces have to be reindexed.
    Inputs:
        *input_surf: surface mesh to be transformed.
        *input_orig: freesurfer orig.mgz.
        *input_deform: deformation (coordinate mapping).
        *input_target: target volume.
        *hemi: hemisphere.
        *path_output: path where to save output.
        *smooth_iter: number of smoothing iterations applied to final image (if set > 0).
        *sort_faces: get new face array if vertices are cut off during deformation.
        *flip_faces: reverse normal direction of mesh.
        *cleanup: remove intermediate files.
        
    created by Daniel Haenelt
    Date created: 06-02-2019          
    Last modified: 22-01-2020
    """
    import os
    import numpy as np
    import nibabel as nb
    import shutil as sh
    from nibabel.freesurfer.io import write_geometry, read_geometry
    from nibabel.affines import apply_affine
    from nipype.interfaces.freesurfer import SampleToSurface
    from nipype.interfaces.freesurfer import SmoothTessellation
    from lib.io.get_filename import get_filename
    from lib.io.mgh2nii import mgh2nii
    from lib.surface.vox2ras import vox2ras

    # set freesurfer path environment
    os.environ["SUBJECTS_DIR"] = path_output

    # freesurfer subject
    tmp = np.random.randint(0, 10, 5)
    tmp_string = ''.join(str(i) for i in tmp)
    sub = "tmp_"+tmp_string

    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # mimic freesurfer folder structure (with some additional folder for intermediate files)
    path_sub = os.path.join(path_output,sub)
    path_mri = os.path.join(path_sub,"mri")
    path_surf = os.path.join(path_sub,"surf")

    os.makedirs(path_sub)
    os.makedirs(path_mri)
    os.makedirs(path_surf)

    # get file extension of orig
    _, name_orig, ext_orig = get_filename(input_orig)

    # name of surface file
    name_surf = os.path.basename(input_surf)

    # copy orig, cmap and input surface to mimicked freesurfer folders
    sh.copyfile(input_surf, os.path.join(path_surf,hemi+".source"))
    if ext_orig != ".mgz":
        mgh2nii(input_orig, path_mri, "mgz")
        os.rename(os.path.join(path_mri,name_orig+".mgz"),os.path.join(path_mri,"orig.mgz"))
    else:
        sh.copyfile(input_orig, os.path.join(path_mri,"orig.mgz"))

    # read surface geometry
    vtx, fac = read_geometry(input_surf)

    # get affine vox2ras-tkr transformation to target volume
    vox2ras_tkr, _ = vox2ras(input_target)
    
    # divide coordinate mapping into its x, y and z components
    cmap_img = nb.load(input_deform)
    cmap_img.header["dim"][0] = 3
    cmap_img.header["dim"][4] = 1

    # apply vox2ras transformation to coordinate mappings
    cmap_array = cmap_img.get_fdata()
    cmap_array = apply_affine(vox2ras_tkr,cmap_array)

    components = ["x", "y", "z"]
    vtx_new = np.zeros([len(vtx),3])
    for i in range(len(components)):
        temp_array = cmap_array[:,:,:,i]
        temp_img = nb.Nifti1Image(temp_array, cmap_img.affine, cmap_img.header)
        nb.save(temp_img,os.path.join(path_mri,components[i]+"_deform.nii"))

        # mri_vol2surf
        sampler = SampleToSurface()
        sampler.inputs.subject_id = sub
        sampler.inputs.reg_header = True
        sampler.inputs.hemi = hemi
        sampler.inputs.source_file = os.path.join(path_mri,components[i]+"_deform.nii")
        sampler.inputs.surface = "source"
        sampler.inputs.sampling_method = "point"
        sampler.inputs.sampling_range = 0
        sampler.inputs.sampling_units = "mm"
        sampler.inputs.interp_method = "nearest" # or trilinear
        sampler.inputs.out_type = "mgh"
        sampler.inputs.out_file = os.path.join(path_surf,hemi+"."+components[i]+"_sampled.mgh")
        sampler.run()
           
        data_img = nb.load(os.path.join(path_surf,hemi+"."+components[i]+"_sampled.mgh"))
        vtx_new[:,i] = np.squeeze(data_img.get_fdata())
    
    if sort_faces:
        
        # get binary mask of slab
        background_array = np.ones(cmap_img.header["dim"][1:4])
        background_array[cmap_img.get_fdata()[:,:,:,0] == 0] = 0
        background_array[cmap_img.get_fdata()[:,:,:,1] == 0] = 0
        background_array[cmap_img.get_fdata()[:,:,:,2] == 0] = 0
    
        background_img = nb.Nifti1Image(background_array, cmap_img.affine, cmap_img.header)
        nb.save(background_img,os.path.join(path_mri,"background.nii"))
        
        # mri_vol2surf (background)
        sampler = SampleToSurface()
        sampler.inputs.subject_id = sub
        sampler.inputs.reg_header = True
        sampler.inputs.hemi = hemi
        sampler.inputs.source_file = os.path.join(path_mri,"background.nii")
        sampler.inputs.surface = "source"
        sampler.inputs.sampling_method = "point"
        sampler.inputs.sampling_range = 0
        sampler.inputs.sampling_units = "mm"
        sampler.inputs.interp_method = "nearest"
        sampler.inputs.out_type = "mgh"
        sampler.inputs.out_file = os.path.join(path_surf,hemi+".background.mgh")
        sampler.run()
        
        # get new indices
        background_list = nb.load(os.path.join(path_surf,hemi+".background.mgh")).get_fdata()
        background_list = np.squeeze(background_list).astype(int)
        
        # only keep vertex indices within the slab
        ind_keep = np.arange(0,len(vtx[:,0]))
        ind_keep[background_list == 0] = -1
        ind_keep = ind_keep[ind_keep != -1]
    
        # get new vertices
        vtx_new = vtx_new[ind_keep,:]
    
        # get new faces
        fac_keep = np.zeros(len(fac[:,0]))
        fac_keep += np.in1d(fac[:,0],ind_keep)
        fac_keep += np.in1d(fac[:,1],ind_keep)
        fac_keep += np.in1d(fac[:,2],ind_keep)
        fac_temp = fac[fac_keep == 3,:]
        fac_new = fac[fac_keep == 3,:]

        # sort new faces
        c_step = 0
        n_step = [10,20,30,40,50,60,70,80,90,100]
        for i in range(len(ind_keep)):
            temp = np.where(ind_keep[i] == fac_temp)
            fac_new[temp] = i
            
            # print status
            counter = np.floor(i / len(ind_keep) * 100).astype(int)
            if counter == n_step[c_step]:
                print("sort faces: "+str(counter)+" %")
                c_step += 1
            
        # remove singularities (vertices without faces)
        fac_counter = 0
        fac_old = fac_new.copy()
        n_singularity = np.zeros(len(vtx_new))
        c_step = 0
        for i in range(len(vtx_new)):
            row, col = np.where(fac_old == i)
         
            n_singularity[i] = len(row)
            if not n_singularity[i]:    
                fac_temp = fac_new.copy()
                fac_temp[fac_temp >= fac_counter] = -1
                fac_temp[fac_temp != -1] = 0
                fac_new += fac_temp
                fac_counter -= 1
         
            # update face counter
            fac_counter += 1
         
            # print status
            counter = np.floor(i / len(vtx_new) * 100).astype(int)
            if counter == n_step[c_step]:
                print("clean vertices: "+str(counter)+" %")
                c_step += 1
        
        # vertices and indices without singularities
        vtx_new = vtx_new[n_singularity != 0]
        ind_keep = ind_keep[n_singularity != 0]
        
        # save index mapping between original and transformed surface
        np.savetxt(os.path.join(path_output, name_surf+"_ind.txt"), ind_keep, fmt='%d')
    else:
        fac_new = fac
 
    # flip faces
    if flip_faces:
        fac_new = np.flip(fac_new, axis=1)
    
    # write new surface
    write_geometry(os.path.join(path_output, name_surf+"_def"), vtx_new, fac_new)

    # smooth surface
    if smooth_iter:
        smooth = SmoothTessellation()
        smooth.inputs.in_file = os.path.join(path_output, name_surf+"_def")
        smooth.inputs.out_file = os.path.join(path_output, name_surf+"_def_smooth")
        smooth.inputs.smoothing_iterations = smooth_iter
        smooth.inputs.disable_estimates = True
        smooth.run()

    # delete intermediate files
    if cleanup:
        sh.rmtree(path_sub, ignore_errors=True)