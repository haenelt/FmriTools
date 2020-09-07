def mask_epi(file_epi, file_t1, file_mask, niter, sigma, file_reg=""):
    """
    This function masks a mean epi image based on a skullstrip mask of the corresponding anatomy.
    The mask is transformed to native epi space via scanner transformation and rigid registration
    of anatomy and epi if no coordinate mapping is given. Finally, holes in the mask are filled, 
    the mask is dilated and a Gaussian filter is applied. The masked epi is saved in the same 
    folder with the prefix p.
    Inputs:
        *file_epi: input mean epi image.
        *file_t1: input of corresponding skullstripped anatomy.
        *file_mask: input of skullstrip mask of the corresponding anatomy.
        *niter: number of dilation iterations.
        *sigma: gaussian smoothing kernel.
        *file_reg: filename of ana -> epi coordinate mapping.
        
    created by Daniel Haenelt
    Date created: 13-02-2019
    Last modified: 07-09-2020
    """
    import os
    import numpy as np
    import nibabel as nb
    import shutil as sh
    from nipype.interfaces.fsl import FLIRT
    from nipype.interfaces.fsl.preprocess import ApplyXFM
    from scipy.ndimage import binary_fill_holes, gaussian_filter
    from scipy.ndimage.morphology import binary_dilation
    from nighres.registration import apply_coordinate_mappings
    from lib.io.get_filename import get_filename
    from lib.cmap.generate_coordinate_mapping import generate_coordinate_mapping
    from lib.cmap.expand_coordinate_mapping import expand_coordinate_mapping
    from lib.registration.get_scanner_transform import get_scanner_transform

    # get paths and filenames
    path_t1, name_t1, _ = get_filename(file_t1)
    path_epi, name_epi, _ = get_filename(file_epi)

    if file_reg:
        _, _, ext_reg = get_filename(file_reg)
    else:
        ext_reg = '.nii.gz'

    # filenames
    file_cmap = os.path.join(path_t1, "cmap.nii")
    file_cmap_reg = os.path.join(path_t1, "cmap_reg"+ext_reg)
    file_cmap_flirt = os.path.join(path_t1, "cmap_flirt.nii.gz")
    file_cmap_def = os.path.join(path_t1, "cmap_def.nii.gz")
    file_ana_reg = os.path.join(path_t1, "ana_reg.nii.gz")
    file_ana_flirt = os.path.join(path_t1, "ana_flirt.nii.gz")
    file_ana_def = os.path.join(path_t1, "ana_def.nii.gz")
    file_mask_def = os.path.join(path_t1, "mask_def.nii.gz")
    file_mask_def2 = os.path.join(path_t1, "mask_def2.nii.gz")
    file_flirt = os.path.join(path_t1, "flirt_matrix.txt")            

    # create new cmap
    cmap = generate_coordinate_mapping(file_epi,
                                       pad=0, 
                                       path_output=None, 
                                       suffix=None, 
                                       time=False, 
                                       write_output=False)
    nb.save(cmap, file_cmap)

    # get initial ana -> epi transformation from existing cmap or header
    if file_reg:
        sh.copyfile(file_reg, file_cmap_reg)        
    else:
        get_scanner_transform(file_t1, file_epi, path_t1, True)
        os.rename(os.path.join(path_t1, name_t1+"_2_"+name_epi+"_scanner.nii.gz"),
                  file_cmap_reg)
    
    # scanner transform peeled t1 to epi
    res = apply_coordinate_mappings(file_t1, # input 
                                    file_cmap_reg, # cmap
                                    interpolation="linear", # nearest or linear
                                    padding="zero", # closest, zero or max
                                    save_data=False,
                                    overwrite=False,
                                    output_dir=None,
                                    file_name=None
                                    )
    nb.save(res["result"], file_ana_reg)

    # flirt t1 to epi
    os.chdir(path_t1)
    flirt = FLIRT()
    flirt.inputs.cost_func = 'corratio'
    flirt.inputs.dof = 6
    flirt.inputs.interp = "trilinear" # trilinear, nearestneighbour, sinc or spline
    flirt.inputs.in_file = file_ana_reg
    flirt.inputs.reference = file_epi
    flirt.inputs.output_type = "NIFTI_GZ"
    flirt.inputs.out_file = file_ana_flirt
    flirt.inputs.out_matrix_file = file_flirt
    flirt.run()

    # apply flirt to generated cmap
    applyxfm = ApplyXFM()
    applyxfm.inputs.in_file = file_cmap
    applyxfm.inputs.reference = file_epi
    applyxfm.inputs.in_matrix_file = file_flirt
    applyxfm.inputs.interp = "trilinear"
    applyxfm.inputs.padding_size = 0
    applyxfm.inputs.output_type = "NIFTI_GZ"
    applyxfm.inputs.out_file = file_cmap_flirt
    applyxfm.inputs.apply_xfm = True
    applyxfm.run() 
    
    # remove outliers and expand
    cmap = nb.load(file_cmap_flirt)
    arr_cmap = cmap.get_fdata()
    
    pts_cmap0 = arr_cmap[0,0,0,0]
    pts_cmap1 = arr_cmap[0,0,0,1]
    pts_cmap2 = arr_cmap[0,0,0,2]

    arr_cmap[arr_cmap == 0] = 0
    arr_cmap[arr_cmap == pts_cmap0] = 0
    arr_cmap[arr_cmap == pts_cmap1] = 0
    arr_cmap[arr_cmap == pts_cmap2] = 0

    output = nb.Nifti1Image(arr_cmap, cmap.affine, cmap.header)
    nb.save(output, file_cmap_flirt)
    
    expand_coordinate_mapping(cmap_in=file_cmap_flirt, 
                              path_output=path_t1,
                              name_output="cmap_flirt", 
                              write_output=True)
        
    # apply flirt cmap to header transformation
    res = apply_coordinate_mappings(file_cmap_reg, # input
                                    file_cmap_flirt,
                                    interpolation="linear", # nearest or linear
                                    padding="zero", # closest, zero or max
                                    save_data=False,
                                    overwrite=False,
                                    output_dir=None,
                                    file_name=None
                                    )
    nb.save(res["result"], file_cmap_def)
    
    # remove outliers and expand
    cmap = nb.load(file_cmap_def)
    arr_cmap = cmap.get_fdata()
    
    pts_cmap0 = arr_cmap[0,0,0,0]
    pts_cmap1 = arr_cmap[0,0,0,1]
    pts_cmap2 = arr_cmap[0,0,0,2]

    arr_cmap[arr_cmap == 0] = 0
    arr_cmap[arr_cmap == pts_cmap0] = 0
    arr_cmap[arr_cmap == pts_cmap1] = 0
    arr_cmap[arr_cmap == pts_cmap2] = 0

    output = nb.Nifti1Image(arr_cmap, cmap.affine, cmap.header)
    nb.save(output, file_cmap_def)
    
    expand_coordinate_mapping(cmap_in=file_cmap_def, 
                              path_output=path_t1,
                              name_output="cmap_def", 
                              write_output=True)
    
    # apply final cmap to t1 and mask
    res = apply_coordinate_mappings(file_t1, # input 
                                    file_cmap_def,
                                    interpolation="linear", # nearest or linear
                                    padding="zero", # closest, zero or max
                                    save_data=False,
                                    overwrite=False,
                                    output_dir=None,
                                    file_name=None
                                    )
    nb.save(res["result"], file_ana_def)
    
    res = apply_coordinate_mappings(file_mask, # input 
                                    file_cmap_def,
                                    interpolation="nearest", # nearest or linear
                                    padding="zero", # closest, zero or max
                                    save_data=False,
                                    overwrite=False,
                                    output_dir=None,
                                    file_name=None
                                    )
    nb.save(res["result"], file_mask_def)
    
    # finalise mask
    mask_img = nb.load(file_mask_def)
    mask_array = mask_img.get_fdata()
    mask_array = binary_fill_holes(mask_array).astype(int) # fill holes in mask
    mask_array = binary_dilation(mask_array, iterations=niter).astype(np.float) # dilate mask
    mask_array = gaussian_filter(mask_array, sigma=sigma) # apply gaussian filter
    
    # write final epi mask
    out_img = nb.Nifti1Image(mask_array, mask_img.affine, mask_img.header)
    nb.save(out_img, file_mask_def2)

    # multiply epi and binary mask
    epi_img = nb.load(file_epi)
    epi_array = epi_img.get_fdata()
    epi_array = epi_array * mask_array # multiply epi and mask
    
    # write masked epi
    out_img = nb.Nifti1Image(epi_array, epi_img.affine, epi_img.header)
    nb.save(out_img, os.path.join(path_epi,"p"+name_epi+".nii"))
