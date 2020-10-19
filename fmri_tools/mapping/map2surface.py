# -*- coding: utf-8 -*-

# python standard library inputs
import os
import shutil as sh

# external inputs
import numpy as np
import nibabel as nb
from nibabel.freesurfer.io import read_geometry
from nipype.interfaces.freesurfer import SampleToSurface
from nipype.interfaces.freesurfer.preprocess import MRIConvert
    

def map2surface(input_surf, input_vol, hemi, path_output, 
                interp_method="nearest", input_white=None, input_ind=None, 
                cleanup=True):
    """ Map to surface

    This function samples data from the input volume to the input surface and 
    optionally maps those values to a target surface if an index file is given.      

    Parameters
    ----------
    input_surf : str
        Surface mesh onto which volume data is sampled.
    input_vol : str
        Volume from which data is sampled.
    hemi : str
        Hemisphere.
    path_output : str
        Path where to save output.
    interp_method : str, optional
        Interpolation method (nearest or trilinear). The default is "nearest".
    input_white : str, optional
        White surface in target surface space (only necessary if index file is 
                                               given). The default is None.
    input_ind : str, optional
        Textfile with mapping of vertex indices to target space. The default is 
        None.
    cleanup : bool, optional
        Remove intermediate files. The default is True.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 06-02-2019      
    Last modified: 19-10-2020

    """
    
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

    # copy input volume as orig.mgz to mimicked freesurfer folder
    if os.path.splitext(os.path.basename(input_vol))[1] != ".mgz":
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
    if os.path.splitext(os.path.basename(input_vol))[1] == ".gz":
        name_vol = os.path.splitext(os.path.splitext(os.path.basename(input_vol))[0])[0]
    else:
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
    sampler.inputs.interp_method = interp_method # nearest or trilinear
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
        nb.save(res_img,os.path.join(path_output,hemi+"."+name_vol+"_"+name_surf+"_def_trans.mgh"))
    else:
        # write sampled data in epi space
        sh.copyfile(os.path.join(path_surf,hemi+"."+"sampled.mgh"),
                    os.path.join(path_output,hemi+"."+name_vol+"_"+name_surf+"_def.mgh"))
    
    # delete intermediate files
    if cleanup:
        sh.rmtree(path_sub, ignore_errors=True)
