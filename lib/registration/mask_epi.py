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
    Last modified: 28-02-2019
    """
    import os
    import numpy as np
    import nibabel as nb
    from nipype.interfaces.fsl import FLIRT
    from nipype.interfaces.fsl.preprocess import ApplyXFM
    from scipy.ndimage import binary_fill_holes, gaussian_filter
    from scipy.ndimage.morphology import binary_dilation
    from nighres.registration import apply_coordinate_mappings
    from lib.registration.get_scanner_transform import get_scanner_transform
    
    # get paths
    path_t1 = os.path.dirname(t1)
    path_epi = os.path.dirname(epi)

    # get filenames
    if os.path.splitext(os.path.basename(t1))[1] == '.gz':
        name_t1 = os.path.splitext(os.path.splitext(os.path.basename(t1))[0])[0]
    else:
        name_t1 = os.path.splitext(os.path.basename(t1))[0]

    if os.path.splitext(os.path.basename(mask))[1] == '.gz':
        name_mask = os.path.splitext(os.path.splitext(os.path.basename(mask))[0])[0]
    else:
        name_mask = os.path.splitext(os.path.basename(mask))[0]

    if os.path.splitext(os.path.basename(epi))[1] == '.gz':
        name_epi = os.path.splitext(os.path.splitext(os.path.basename(epi))[0])[0]
    else:
        name_epi = os.path.splitext(os.path.basename(epi))[0]

    # get scanner transform
    get_scanner_transform(t1, epi, path_t1, False)

    # scanner transform peeled t1 to epi
    apply_coordinate_mappings(t1, # input 
                              os.path.join(path_t1,name_t1+"_2_"+name_epi+"_scanner.nii"), # cmap
                              interpolation = "linear", # nearest or linear
                              padding = "zero", # closest, zero or max
                              save_data = True, # save output data to file (boolean)
                              overwrite = True, # overwrite existing results (boolean)
                              output_dir = path_t1, # output directory
                              file_name = name_t1 # base name with file extension for output
                              )

    # scanner transform mask to epi
    apply_coordinate_mappings(mask, # input 
                              os.path.join(path_t1,name_t1+"_2_"+name_epi+"_scanner.nii"), # cmap
                              interpolation = "linear", # nearest or linear
                              padding = "zero", # closest, zero or max
                              save_data = True, # save output data to file (boolean)
                              overwrite = True, # overwrite existing results (boolean)
                              output_dir = path_t1, # output directory
                              file_name = name_mask # base name with file extension for output
                              )

    # flirt t1 to epi
    os.chdir(path_t1)
    flirt = FLIRT()
    flirt.inputs.cost_func = "mutualinfo"
    flirt.inputs.dof = 6
    flirt.inputs.interp = "trilinear" # trilinear, nearestneighbour, sinc or spline
    flirt.inputs.in_file = os.path.join(path_t1,name_t1+"_def-img.nii.gz")
    flirt.inputs.reference = epi
    flirt.inputs.output_type = "NIFTI_GZ"
    flirt.inputs.out_file = os.path.join(path_t1,name_t1+"_def-img2.nii.gz")
    flirt.inputs.out_matrix_file = os.path.join(path_t1,"flirt_matrix.txt")
    flirt.run()

    # flirt mask to epi
    applyxfm = ApplyXFM()
    applyxfm.inputs.in_file = os.path.join(path_t1,name_mask+"_def-img.nii.gz")
    applyxfm.inputs.reference = epi
    applyxfm.inputs.in_matrix_file = os.path.join(path_t1,"flirt_matrix.txt")
    applyxfm.inputs.interp = "trilinear"
    applyxfm.inputs.padding_size = 0
    applyxfm.inputs.output_type = "NIFTI_GZ"
    applyxfm.inputs.out_file = os.path.join(path_t1,name_mask+"_def-img2.nii.gz")
    applyxfm.inputs.apply_xfm = True
    applyxfm.run() 

    # finalise mask
    mask_img = nb.load(os.path.join(path_t1,name_mask+"_def-img2.nii.gz"))
    mask_array = mask_img.get_fdata()
    mask_array = binary_fill_holes(mask_array).astype(int) # fill holes in mask
    mask_array = binary_dilation(mask_array, iterations=niter).astype(np.float) # dilate mask
    mask_array = gaussian_filter(mask_array, sigma=sigma) # apply gaussian filter
    
    # write final epi mask
    out_img = nb.Nifti1Image(mask_array, mask_img.affine, mask_img.header)
    nb.save(out_img,os.path.join(path_t1,name_mask+"_def-img3.nii.gz"))

    # multiply epi and binary mask
    epi_img = nb.load(epi)
    epi_array = epi_img.get_fdata()
    epi_array = epi_array * mask_array # multiply epi and mask
    
    # write masked epi
    out_img = nb.Nifti1Image(epi_array, epi_img.affine, epi_img.header)
    nb.save(out_img,os.path.join(path_epi,"p"+name_epi+".nii"))
