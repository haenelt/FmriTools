def map2surface(input_surf, input_vol, hemi, path_output, input_white=None, input_ind=None, cleanup=True):
    """
    This function samples data from the input volume to the input surface and optionally maps those
    values to a target surface if an index file is given.
    Inputs:
        *input_surf: surface mesh onto which volume data is sampled.
        *input_vol: volume from which data is sampled.
        *hemi: hemisphere.
        *input_white: white surface in target surface space (only necessary if index file is given).
        *input_ind: textfile with mapping of vertex indices to target space.
        *path_output: path where to save output.
        *cleanup: remove intermediate files.
            
    created by Daniel Haenelt
    Date created: 06-02-2019      
    Last modified: 18-02-2019
    """
    import os
    import numpy as np
    import nibabel as nb
    import shutil as sh
    from nibabel.freesurfer.io import read_geometry
    from nipype.interfaces.freesurfer import SampleToSurface
    from nipype.interfaces.freesurfer.preprocess import MRIConvert

    # set freesurfer path environment
    os.environ["SUBJECTS_DIR"] = path_output

    # freesurfer subject
    tmp = np.random.randint(0, 10, 5)
    tmp_string = ''.join(str(i) for i in tmp)
    sub = "tmp_"+tmp_string

    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # mimic freesurfer folder structure (with some additional folder for intermediate files)
    path_sub = os.path.join(path_output,sub)
    path_mri = os.path.join(path_sub,"mri")
    path_surf = os.path.join(path_sub,"surf")

    os.makedirs(path_sub)
    os.makedirs(path_mri)
    os.makedirs(path_surf)

    # copy input volume as orig.mgz to mimicked freesurfer folder
    if os.path.splitext(os.path.basename(input_vol))[1] is not ".mgz":
        mc = MRIConvert()
        mc.inputs.in_file = input_vol
        mc.inputs.out_file = os.path.join(path_mri,"orig.mgz")
        mc.inputs.out_type = 'mgz'
        mc.run()
    else:
        sh.copyfile(input_vol,os.path.join(path_mri,"orig.mgz"))

    # copy input surface to mimicked freesurfer folder
    sh.copyfile(input_surf, os.path.join(path_surf,hemi+".source"))

    # input volume file name
    name_vol = os.path.splitext(os.path.basename(input_vol))[0]
    name_surf = os.path.basename(input_surf).split('.')[1]

    # mri_vol2surf
    sampler = SampleToSurface()
    sampler.inputs.subject_id = sub
    sampler.inputs.reg_header = True
    sampler.inputs.hemi = hemi
    sampler.inputs.source_file = input_vol
    sampler.inputs.surface = "source"
    sampler.inputs.sampling_method = "point"
    sampler.inputs.sampling_range = 0
    sampler.inputs.sampling_units = "mm"
    sampler.inputs.interp_method = "nearest" # or trilinear
    sampler.inputs.out_type = "mgh"
    sampler.inputs.out_file = os.path.join(path_surf,hemi+"."+"sampled.mgh")
    sampler.run()

    if input_ind:
        # read ind
        ind_orig = np.loadtxt(input_ind, dtype=int)
    
        # read white
        vtx_orig, _ = read_geometry(input_white)
    
        # read sampled morph data
        vals_img = nb.load(os.path.join(path_surf,hemi+"."+"sampled.mgh"))
        vals_array = vals_img.get_fdata()
    
        # get anzahl der vertices als Nullen in white
        vals_orig = np.zeros([len(vtx_orig[:,0]),1,1])
    
        # setze sampled data da rein
        vals_orig[ind_orig] = vals_array
    
        # write sampled data in anatomical space      
        vals_img.header["dims"][0] = len(vals_orig[:,0])
        vals_img.header["Mdc"] = np.eye(3)
        res_img = nb.Nifti1Image(vals_orig,vals_img.affine,vals_img.header)
        nb.save(res_img,os.path.join(path_output,hemi+"."+name_vol+"_"+name_surf+"_def_ana.mgh"))
    else:
        # write sampled data in epi space
        sh.copyfile(os.path.join(path_surf,hemi+"."+"sampled.mgh"),
                    os.path.join(path_output,hemi+"."+name_vol+"_"+name_surf+"_def_epi.mgh"))
    
    # delete intermediate files
    if cleanup:
        sh.rmtree(path_sub, ignore_errors=True)
