def mask_epi(epi, t1, mask, niter, sigma):
    """
    This function masks a mean epi image based on a skullstrip mask of the corresponding anatomy.
    The mask is transformed to native epi space via scanner transformation and rigid registration
    of anatomy and epi. Finally, holes in the mask are filled, the mask is dilated and a Gaussian
    filter is applied. The masked epi is saved in the same folder with the prefix p.
    Inputs:
        *epi: input mean epi image.
        *t1: input of corresponding skullstripped anatomy.
        *mask: input of skullstrip mask of the corresponding anatomy.
        *niter: number of dilation iterations.
        *sigma: gaussian smoothing kernel.
        
    created by Daniel Haenelt
    Date created: 13-02-2019
    Last modified: 26-06-2020
    """
    import os
    import numpy as np
    import nibabel as nb
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
    path_t1, name_t1, _ = get_filename(t1)
    path_epi, name_epi, _ = get_filename(epi)
    _, name_mask, _ = get_filename(mask)

    # create new cmap
    cmap = generate_coordinate_mapping(epi, 
                                       pad=0, 
                                       path_output=None, 
                                       suffix=None, 
                                       time=False, 
                                       write_output=False)
    nb.save(cmap, os.path.join(path_t1, "cmap.nii"))

    # get scanner transform
    get_scanner_transform(t1, epi, path_t1, False)

    # scanner transform peeled t1 to epi
    res = apply_coordinate_mappings(t1, # input 
                                    os.path.join(path_t1,name_t1+"_2_"+name_epi+"_scanner.nii"), # cmap
                                    interpolation="linear", # nearest or linear
                                    padding="zero", # closest, zero or max
                                    save_data=False,
                                    overwrite=False,
                                    output_dir=None,
                                    file_name=None
                                    )
    nb.save(res["result"], os.path.join(path_t1, name_t1+"_header.nii.gz"))

    # flirt t1 to epi
    os.chdir(path_t1)
    flirt = FLIRT()
    flirt.inputs.cost_func = "corratio"
    flirt.inputs.dof = 6
    flirt.inputs.interp = "trilinear" # trilinear, nearestneighbour, sinc or spline
    flirt.inputs.in_file = os.path.join(path_t1, name_t1+"_header.nii.gz")
    flirt.inputs.reference = epi
    flirt.inputs.output_type = "NIFTI_GZ"
    flirt.inputs.out_file = os.path.join(path_t1, name_t1+"_header_flirt.nii.gz")
    flirt.inputs.out_matrix_file = os.path.join(path_t1,"flirt_matrix.txt")
    flirt.run()

    # apply flirt to generated cmap
    applyxfm = ApplyXFM()
    applyxfm.inputs.in_file = os.path.join(path_t1,"cmap.nii")
    applyxfm.inputs.reference = epi
    applyxfm.inputs.in_matrix_file = os.path.join(path_t1,"flirt_matrix.txt")
    applyxfm.inputs.interp = "trilinear"
    applyxfm.inputs.padding_size = 0
    applyxfm.inputs.output_type = "NIFTI_GZ"
    applyxfm.inputs.out_file = os.path.join(path_t1,"cmap_flirt.nii.gz")
    applyxfm.inputs.apply_xfm = True
    applyxfm.run() 
    
    # remove outliers and expand
    cmap = nb.load(os.path.join(path_t1, "cmap_flirt.nii.gz"))
    arr_cmap = cmap.get_fdata()
    
    pts_cmap0 = arr_cmap[0,0,0,0]
    pts_cmap1 = arr_cmap[0,0,0,1]
    pts_cmap2 = arr_cmap[0,0,0,2]

    arr_cmap[arr_cmap == 0] = 0
    arr_cmap[arr_cmap == pts_cmap0] = 0
    arr_cmap[arr_cmap == pts_cmap1] = 0
    arr_cmap[arr_cmap == pts_cmap2] = 0

    output = nb.Nifti1Image(arr_cmap, cmap.affine, cmap.header)
    nb.save(output, os.path.join(path_t1, "cmap_flirt.nii.gz"))
    
    expand_coordinate_mapping(cmap_in=os.path.join(path_t1,"cmap_flirt.nii.gz"), 
                              path_output=path_t1,
                              name_output="cmap_flirt", 
                              write_output=True)
        
    # apply flirt cmap to header transformation
    res = apply_coordinate_mappings(os.path.join(path_t1,name_t1+"_2_"+name_epi+"_scanner.nii"), # input 
                                    os.path.join(path_t1,"cmap_flirt.nii.gz"),
                                    interpolation="linear", # nearest or linear
                                    padding="zero", # closest, zero or max
                                    save_data=False,
                                    overwrite=False,
                                    output_dir=None,
                                    file_name=None
                                    )
    nb.save(res["result"], os.path.join(path_t1, "cmap_def.nii.gz"))
    
    # remove outliers and expand
    cmap = nb.load(os.path.join(path_t1, "cmap_def.nii.gz"))
    arr_cmap = cmap.get_fdata()
    
    pts_cmap0 = arr_cmap[0,0,0,0]
    pts_cmap1 = arr_cmap[0,0,0,1]
    pts_cmap2 = arr_cmap[0,0,0,2]

    arr_cmap[arr_cmap == 0] = 0
    arr_cmap[arr_cmap == pts_cmap0] = 0
    arr_cmap[arr_cmap == pts_cmap1] = 0
    arr_cmap[arr_cmap == pts_cmap2] = 0

    output = nb.Nifti1Image(arr_cmap, cmap.affine, cmap.header)
    nb.save(output, os.path.join(path_t1, "cmap_def.nii.gz"))
    
    expand_coordinate_mapping(cmap_in=os.path.join(path_t1,"cmap_def.nii.gz"), 
                              path_output=path_t1,
                              name_output="cmap_def", 
                              write_output=True)
    
    # apply final cmap to t1 and mask
    res = apply_coordinate_mappings(t1, # input 
                                    os.path.join(path_t1,"cmap_def.nii.gz"),
                                    interpolation="linear", # nearest or linear
                                    padding="zero", # closest, zero or max
                                    save_data=False,
                                    overwrite=False,
                                    output_dir=None,
                                    file_name=None
                                    )
    nb.save(res["result"], os.path.join(path_t1,name_t1+"_def.nii"))
    
    res = apply_coordinate_mappings(mask, # input 
                                    os.path.join(path_t1,"cmap_def.nii.gz"),
                                    interpolation="nearest", # nearest or linear
                                    padding="zero", # closest, zero or max
                                    save_data=False,
                                    overwrite=False,
                                    output_dir=None,
                                    file_name=None
                                    )
    nb.save(res["result"], os.path.join(path_t1,"mask_def.nii"))
    
    # finalise mask
    mask_img = nb.load(os.path.join(path_t1,"mask_def.nii"))
    mask_array = mask_img.get_fdata()
    mask_array = binary_fill_holes(mask_array).astype(int) # fill holes in mask
    mask_array = binary_dilation(mask_array, iterations=niter).astype(np.float) # dilate mask
    mask_array = gaussian_filter(mask_array, sigma=sigma) # apply gaussian filter
    
    # write final epi mask
    out_img = nb.Nifti1Image(mask_array, mask_img.affine, mask_img.header)
    nb.save(out_img,os.path.join(path_t1, "mask_def2.nii"))

    # multiply epi and binary mask
    epi_img = nb.load(epi)
    epi_array = epi_img.get_fdata()
    epi_array = epi_array * mask_array # multiply epi and mask
    
    # write masked epi
    out_img = nb.Nifti1Image(epi_array, epi_img.affine, epi_img.header)
    nb.save(out_img,os.path.join(path_epi,"p"+name_epi+".nii"))