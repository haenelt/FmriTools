def deform_surface(input_surf, input_orig, input_deform, input_target, hemi, path_output, cleanup=True):
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
        *cleanup: remove intermediate files.
        
    created by Daniel Haenelt
    Date created: 06-02-2019          
    Last modified: 07-02-2019
    """
    import os
    import subprocess
    import numpy as np
    import nibabel as nb
    import shutil as sh
    from nibabel.freesurfer.io import write_geometry, read_geometry
    from nibabel.affines import apply_affine
    from nipype.interfaces.freesurfer import SampleToSurface
    from nipype.interfaces.freesurfer import SmoothTessellation

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

    # copy orig and input surface to mimicked freesurfer folders
    sh.copyfile(input_surf, os.path.join(path_surf,hemi+".source"))
    sh.copyfile(input_orig, os.path.join(path_mri,"orig.mgz"))

    # read surface geometry
    vtx, fac = read_geometry(input_surf)

    # get affine vox2ras-tkr transformation to target volume
    transformation = subprocess.check_output(['mri_info', input_target, '--{}'.format("vox2ras-tkr")]).decode()
    num_transformation = [[float(x) for x in line.split()] for line in transformation.split('\n') if len(line)>0]
    vox2ras_tkr = np.array(num_transformation)
    
    # divide coordinate mapping into its x, y and z components
    cmap_img = nb.load(input_deform)
    cmap_img.header["dim"][0] = 3
    cmap_img.header["dim"][4] = 1
    cmap_img.header["pixdim"][3] = 1

    # apply vox2ras transformation to coordinate mappings
    cmap_array = cmap_img.get_fdata()
    cmap_array = apply_affine(vox2ras_tkr,cmap_array)
    cmap_background = apply_affine(vox2ras_tkr,[0,0,0])
    cmap_array[:,:,:,0][cmap_array[:,:,:,0]==cmap_background[0]] = 0
    cmap_array[:,:,:,1][cmap_array[:,:,:,1]==cmap_background[1]] = 0
    cmap_array[:,:,:,2][cmap_array[:,:,:,2]==cmap_background[2]] = 0

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
    
    # get new indices
    ind_keep = np.arange(0,len(vtx[:,0]))
    ind_keep = ind_keep[np.sum(vtx_new, axis=1) != 0]

    # get new vertices
    vtx_new = vtx_new[ind_keep]

    # get new faces
    fac_keep = np.zeros(len(fac[:,0]))
    fac_keep += np.in1d(fac[:,0],ind_keep)
    fac_keep += np.in1d(fac[:,1],ind_keep)
    fac_keep += np.in1d(fac[:,2],ind_keep)
    fac_temp = fac[fac_keep == 3,:]
    fac_new = fac[fac_keep == 3,:]

    # sort new faces
    for i in range(len(ind_keep)):
        temp = np.where(ind_keep[i] == fac_temp)
        fac_new[temp] = i

    # write new surface
    write_geometry(os.path.join(path_surf,hemi+".transformed"), vtx_new, fac_new)

    # smooth surface
    smooth = SmoothTessellation()
    smooth.inputs.in_file = os.path.join(path_surf,hemi+".transformed")
    smooth.inputs.out_file = os.path.join(path_output,os.path.basename(input_surf)+"_def")
    smooth.inputs.smoothing_iterations = 10
    smooth.inputs.disable_estimates = True
    smooth.run()

    # save index mapping between original and transformed surface
    np.savetxt(os.path.join(path_output,os.path.basename(input_surf)+"_ind.txt"), ind_keep, fmt='%d')

    # delete intermediate files
    if cleanup:
        sh.rmtree(path_sub, ignore_errors=True)