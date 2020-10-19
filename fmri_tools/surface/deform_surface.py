# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys
import itertools
import shutil as sh

# external inputs
import numpy as np
import nibabel as nb
from nibabel.freesurfer.io import write_geometry, read_geometry
from nibabel.affines import apply_affine
from nipype.interfaces.freesurfer import SampleToSurface
from nipype.interfaces.freesurfer import SmoothTessellation
from gbb.utils.vox2ras import vox2ras

# local inputs
from fmri_tools.io.get_filename import get_filename
from fmri_tools.io.mgh2nii import mgh2nii


def deform_surface(input_surf, input_orig, input_deform, input_target, 
                   path_output, input_mask=None, interp_method="nearest", 
                   smooth_iter=0, flip_faces=False, cleanup=True):
    """ Deform surface

    This function deforms a surface mesh in freesurfer convention using a 
    coordinate map containing voxel coordinates. The computation takes quite a 
    while because in the case of removed vertices, i.e. if a mask is given as 
    input, the remaining faces are reindexed.    

    Parameters
    ----------
    input_surf : str
        Surface mesh to be transformed.
    input_orig : str
        Freesurfer orig.mgz.
    input_deform : str
        Deformation (coordinate mapping).
    input_target : str
        Target volume.
    path_output : str
        Path where to save output.
    input_mask : str, optional
        Mask volume. The default is None.
    interp_method : str, optional
        Interpolation method (nearest or trilinear). The default is "nearest".
    smooth_iter : int, optional
        Number of smoothing iterations applied to final image (if set > 0). The 
        default is 0.
    flip_faces : bool, optional
        Reverse normal direction of mesh. The default is False.
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

    # mimic freesurfer folder structure (with some additional folder for 
    # intermediate files)
    path_sub = os.path.join(path_output,sub)
    path_mri = os.path.join(path_sub,"mri")
    path_surf = os.path.join(path_sub,"surf")

    os.makedirs(path_sub)
    os.makedirs(path_mri)
    os.makedirs(path_surf)

    # get filenames
    _, name_orig, ext_orig = get_filename(input_orig)
    _, hemi, name_surf = get_filename(input_surf)
    name_surf = name_surf.replace(".","")
    
    # check filename
    if not hemi == "lh" and not hemi == "rh":
        sys.exit("Could not identify hemi from filename!")

    # copy orig, cmap and input surface to mimicked freesurfer folders
    sh.copyfile(input_surf, os.path.join(path_surf,hemi+".source"))
    if ext_orig != ".mgz":
        mgh2nii(input_orig, path_mri, "mgz")
        os.rename(os.path.join(path_mri,name_orig+".mgz"),
                  os.path.join(path_mri,"orig.mgz"))
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
        file_temp = os.path.join(path_mri,components[i]+"_deform.nii")
        file_sampled = os.path.join(path_surf,hemi+"."+components[i]+"_sampled.mgh")
        
        # get target volume
        temp_array = cmap_array[:,:,:,i]
        temp_img = nb.Nifti1Image(temp_array, cmap_img.affine, cmap_img.header)
        nb.save(temp_img, file_temp)

        # mri_vol2surf
        sampler = SampleToSurface()
        sampler.inputs.subject_id = sub
        sampler.inputs.reg_header = True
        sampler.inputs.hemi = hemi
        sampler.inputs.source_file = file_temp
        sampler.inputs.surface = "source"
        sampler.inputs.sampling_method = "point"
        sampler.inputs.sampling_range = 0
        sampler.inputs.sampling_units = "mm"
        sampler.inputs.interp_method = interp_method
        sampler.inputs.out_type = "mgh"
        sampler.inputs.out_file = file_sampled
        sampler.run()
           
        data_img = nb.load(file_sampled)
        vtx_new[:,i] = np.squeeze(data_img.get_fdata())
    
    if input_mask:
        file_background = os.path.join(path_surf,hemi+".background.mgh")
                
        # mri_vol2surf (background)
        sampler = SampleToSurface()
        sampler.inputs.subject_id = sub
        sampler.inputs.reg_header = True
        sampler.inputs.hemi = hemi
        sampler.inputs.source_file = input_mask
        sampler.inputs.surface = "source"
        sampler.inputs.sampling_method = "point"
        sampler.inputs.sampling_range = 0
        sampler.inputs.sampling_units = "mm"
        sampler.inputs.interp_method = "nearest"
        sampler.inputs.out_type = "mgh"
        sampler.inputs.out_file = file_background
        sampler.run()
        
        # get new indices
        background_list = nb.load(file_background).get_fdata()
        background_list = np.squeeze(background_list).astype(int)
        
        # only keep vertex indices within the slab
        ind_keep = np.arange(len(vtx))
        ind_keep = ind_keep[background_list != 0]

        # get indices which will be removed
        ind_tmp = np.arange(len(vtx))
        ind_remove = list(set(ind_tmp) - set(ind_keep))
        ind_remove = sorted(ind_remove, reverse=True)

        # get new vertices
        vtx_new = vtx_new[ind_keep,:]

        # get new faces
        fac_keep = np.zeros(len(fac))
        fac_keep += np.in1d(fac[:,0], ind_keep)
        fac_keep += np.in1d(fac[:,1], ind_keep)
        fac_keep += np.in1d(fac[:,2], ind_keep)
        fac_new = fac[fac_keep == 3,:]

        # reindex faces
        loop_status = 0
        loop_length = len(ind_remove)
        for i in range(loop_length):
                
            # print status
            counter = np.floor(i / loop_length * 100)
            if counter != loop_status:
                print("sort faces: "+str(counter)+" %")
                loop_status = counter
                
            tmp = fac_new[fac_new >= ind_remove[i]] - 1
            fac_new[fac_new >= ind_remove[i]] = tmp

        # get indices which will be cleaned
        ind_vtx_new = np.arange(len(vtx_new))
        ind_fac_new = list(itertools.chain(*fac_new))
        ind_fac_new = list(set(ind_fac_new))
        ind_remove = list(set(ind_vtx_new) - set(ind_fac_new))
        ind_remove = sorted(ind_remove, reverse=True)
        
        # remove singularities (vertices without faces)
        loop_status = 0
        loop_length = len(ind_remove)
        for i in range(loop_length):
            
            # print status
            counter = np.floor(i / loop_length * 100)
            if counter != loop_status:
                print("clean faces: "+str(counter)+" %")
                loop_status = counter
                       
            # remove vertex and index
            vtx_new = np.delete(vtx_new, ind_remove[i], 0)
            ind_keep = np.delete(ind_keep, ind_remove[i], 0)                

            # sort faces
            tmp = fac_new[fac_new >= ind_remove[i]] - 1
            fac_new[fac_new >= ind_remove[i]] = tmp
        
        # save index mapping between original and transformed surface
        np.savetxt(os.path.join(path_output, hemi+"."+name_surf+"_ind.txt"), ind_keep, fmt='%d')
    else:
        fac_new = fac
 
    # flip faces
    if flip_faces:
        fac_new = np.flip(fac_new, axis=1)
    
    # write new surface
    write_geometry(os.path.join(path_output, hemi+"."+name_surf+"_def"), 
                   vtx_new, fac_new)

    # smooth surface
    if smooth_iter:
        smooth = SmoothTessellation()
        smooth.inputs.in_file = os.path.join(path_output, hemi+"."+name_surf+"_def")
        smooth.inputs.out_file = os.path.join(path_output, hemi+"."+name_surf+"_def_smooth")
        smooth.inputs.smoothing_iterations = smooth_iter
        smooth.inputs.disable_estimates = True
        smooth.run()

    # delete intermediate files
    if cleanup:
        sh.rmtree(path_sub, ignore_errors=True)
