def calculate_equivolumetric_epi2(input_white, input_pial, input_vol, path_output, n_layers, 
                                  r=[0.4,0.4,0.4], n_iter=2):
    """
    This function computes equivolumetric layers in volume space from input pial and white surfaces
    in freesurfer format. The input surfaces do not have to cover the whole brain. Number of 
    vertices and indices do not have to correspond between surfaces.
    Inputs:
        *input_white: filename of white surface.
        *input_pial: filename of pial surface.
        *input_vol: filename of reference volume.
        *path_output: path where output is written.
        *n_layers: number of generated layers + 1.
        *r: array of new voxel sizes for reference volume upsampling.
        *n_iter: number of surface upsampling iterations.
    
    created by Daniel Haenelt
    Date created: 17-12-2019
    Last modified: 27-05-2020
    """
    import sys
    import os
    import numpy as np
    import nibabel as nb
    from nibabel.affines import apply_affine
    from nibabel.freesurfer.io import read_geometry
    from skimage import measure
    from nighres.surface import probability_to_levelset
    from nighres.laminar import volumetric_layering
    from scipy.ndimage.morphology import binary_fill_holes
    from collections import Counter
    from lib.utils.upsample_volume import upsample_volume
    from lib.surface.vox2ras import vox2ras
    from lib.surface.upsample_surf_mesh import upsample_surf_mesh
    
    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)
    
    # get hemi from filename
    hemi = os.path.splitext(os.path.basename(input_white))[0]
    if not hemi == "lh" and not hemi == "rh":
        sys.exit("Could not identify hemi from filename!")
    
    # new filenames in output folder
    res_white = os.path.join(path_output,hemi+".white")
    res_pial = os.path.join(path_output,hemi+".pial")
    res_vol = os.path.join(path_output,"epi_upsampled.nii")
    
    # upsample reference volume and input surface
    upsample_volume(input_vol, res_vol, dxyz=r, rmode="Cu")    
    upsample_surf_mesh(input_white, res_white, n_iter, "linear")
    upsample_surf_mesh(input_pial, res_pial, n_iter, "linear")
    
    # get affine ras2vox-tkr transformation to reference volume
    _, ras2vox_tkr = vox2ras(res_vol)
    
    # load surface
    vtx_white, fac_white = read_geometry(res_white) 
    vtx_pial, _ = read_geometry(res_pial)
    
    # load volume
    vol = nb.load(res_vol)
    
    # apply ras2vox to coords
    vtx_white = apply_affine(ras2vox_tkr, vtx_white).astype(int)
    vtx_pial = apply_affine(ras2vox_tkr, vtx_pial).astype(int)
    
    # surfaces to lines in volume
    white_array = np.zeros(vol.header["dim"][1:4])
    white_array[vtx_white[:,0],vtx_white[:,1],vtx_white[:,2]] = 1   
    white = nb.Nifti1Image(white_array, vol.affine, vol.header)
    
    pial_array = np.zeros(vol.header["dim"][1:4])
    pial_array[vtx_pial[:,0],vtx_pial[:,1],vtx_pial[:,2]] = 1
    pial = nb.Nifti1Image(pial_array, vol.affine, vol.header)
    
    """
    make wm
    """
    white_label_array = np.zeros_like(white_array)
    for i in range(np.shape(white_label_array)[2]):
        white_label_array[:,:,i] = binary_fill_holes(white_array[:,:,i])
    white_label_array = white_label_array - white_array
    white_label_array = measure.label(white_label_array, connectivity=1)
    white_label_flatten = np.ndarray.flatten(white_label_array)
    white_label_flatten = white_label_flatten[white_label_flatten > 0]
    label_number = Counter(white_label_flatten).most_common(1)[0][0]
    white_label_array[white_label_array != label_number] = 0
    white_label_array[white_label_array > 0] = 1    
    white_label = nb.Nifti1Image(white_label_array, vol.affine, vol.header)
    
    """
    make csf
    """
    pial_label_array = np.zeros_like(pial_array)
    for i in range(np.shape(pial_label_array)[2]):
        pial_label_array[:,:,i] = binary_fill_holes(pial_array[:,:,i])
    pial_label_array = pial_label_array - pial_array
    pial_label_array = measure.label(pial_label_array, connectivity=1)   
    pial_label_flatten = np.ndarray.flatten(pial_label_array)
    pial_label_flatten = pial_label_flatten[pial_label_flatten > 0]
    label_number = Counter(pial_label_flatten).most_common(1)[0][0]
    pial_label_array[pial_label_array != label_number] = 0
    pial_label_array[pial_label_array > 0] = 1    
    pial_label = nb.Nifti1Image(pial_label_array, vol.affine, vol.header)
       
    """
    make gm
    """
    ribbon_label_array = pial_label_array.copy()
    ribbon_label_array = pial_label_array - white_label_array
    ribbon_label_array[ribbon_label_array != 1] = 0
    ribbon_label = nb.Nifti1Image(ribbon_label_array, vol.affine, vol.header)
    
    """
    layers
    """
    csf_level = probability_to_levelset(pial_label)
    wm_level = probability_to_levelset(white_label)
    
    volumetric_layering(wm_level["result"], 
                        csf_level["result"], 
                        n_layers=n_layers, 
                        topology_lut_dir=None,
                        save_data=True, 
                        overwrite=True, 
                        output_dir=path_output, 
                        file_name="epi")
    
    # write niftis
    nb.save(white, os.path.join(path_output,"wm_line.nii"))
    nb.save(pial, os.path.join(path_output,"csf_line.nii"))
    nb.save(white_label,os.path.join(path_output,"wm_label.nii"))
    nb.save(pial_label,os.path.join(path_output,"csf_label.nii"))
    nb.save(ribbon_label, os.path.join(path_output,"gm_label.nii"))
