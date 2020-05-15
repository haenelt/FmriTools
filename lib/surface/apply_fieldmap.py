def apply_fieldmap(file_fmap_magn, file_fmap_phase, file_epi, file_epi_moco, file_surf,
                   delta_te=1.02, smooth=2.5, udir="y-", bw=16.304, nerode=1, cleanup=True):
    """
    This function computes a deformation field from a fieldmap acquisition and applies the inverse
    transformation to the undistorted surface. The following steps are performed:
        1. get median time series
        2. skullstrip epi
        3. register fieldmap to epi
        4. mask fieldmap
        5. prepare field
        6. get deforamtion field
        7. apply inverse deformation to surfaces.
        8. remove intermediate files (optional).
        
    To run the script, FSL and Freesurfer have to be in the PATH environment. The basenames of the 
    surface files should be in freesurfer convention with the hemisphere indicated as prefix.
    Inputs:
        *fiele_fmap_magn: fieldmap magnitude image.
        *file_fmap_phase: fieldmap phase difference image.
        *file_epi: filename of raw time series.
        *file_epi_moco: filname of motion corrected time series.
        *file_surf: list of surface filnames.
        *delta_te: echo time difference of fieldmap in ms.
        *smooth: smoothing kernel for fieldmap unmasking.
        *udir: direction for fieldmap unmasking.
        *bw: BandwidthPerPixelPhaseEncode in Hz/px.
        *nerode: number of skullstrip mask eroding iterations.
        *cleanup: removes temporary files at the end of the script (boolean).
    
    created by Daniel Haenelt
    Date created: 31-01-2020
    Last modified: 13-05-2020
    """
    import os
    import numpy as np
    import nibabel as nb
    from nipype.interfaces import fsl
    from lib.skullstrip.skullstrip_epi import skullstrip_epi
    from lib.io.get_filename import get_filename
    from lib.cmap.generate_coordinate_mapping import generate_coordinate_mapping
    from lib.surface.deform_surface import deform_surface
    
    # prepare path and filename
    path_fmap0, name_fmap0, ext_fmap0 = get_filename(file_fmap_magn)
    path_fmap1, name_fmap1, ext_fmap1 = get_filename(file_fmap_phase)
    path_data, name_data, ext_data = get_filename(file_epi)
    path_udata, name_udata, ext_udata = get_filename(file_epi_moco)
    
    # filename with file extension
    name_fmap0 += ext_fmap0
    name_fmap1 += ext_fmap1
    name_data += ext_data
    name_udata += ext_udata
    
    # change directory to fieldmap directory
    os.chdir(path_fmap0)
    
    # get matrix size in phase encoding direction from uncorrected epi
    data = nb.load(file_epi)
    phase_encode = data.header.get_dim_info()[1]
    ImageMatrixPhaseEncode = data.header["dim"][phase_encode+1]
    
    # calculate median epi
    udata = nb.load(file_epi_moco)
    arr_udata = udata.get_fdata()
    arr_udata_median = np.median(arr_udata, axis=3)
    udata_median = nb.Nifti1Image(arr_udata_median, udata.affine, udata.header)
    udata_median.header["dim"][0] = 3
    udata_median.header["dim"][4] = 1
    nb.save(udata_median, os.path.join(path_udata, "median_"+name_udata))
    
    # calculate skullstrip mask of that image
    skullstrip_epi(os.path.join(path_udata, "median_"+name_udata),
                   roi_size=10, 
                   scale=0.75, 
                   nerode=1, 
                   ndilate=2, 
                   savemask=True, 
                   cleanup=True)
    
    # erode skullstrip mask
    for j in range(nerode):
        erode = fsl.ErodeImage()
        erode.inputs.in_file = os.path.join(path_udata, "mask_median_"+name_udata)
        erode.inputs.output_type = "NIFTI"
        erode.inputs.out_file = os.path.join(path_udata, "mask_median_"+name_udata)
        erode.run()
    
    # register fmap1 to median epi (fsl.FLIRT)
    flirt = fsl.FLIRT()
    flirt.inputs.cost_func = "mutualinfo"
    flirt.inputs.dof = 6
    flirt.inputs.interp = "trilinear" # trlinear, nearestneighbour, sinc or spline
    flirt.inputs.in_file = file_fmap_magn
    flirt.inputs.reference = os.path.join(path_udata, "median_"+name_udata)
    flirt.inputs.output_type = "NIFTI"
    flirt.inputs.out_file = os.path.join(path_fmap0, "r"+name_fmap0)
    flirt.inputs.out_matrix_file = os.path.join(path_fmap0, "fmap2epi.txt")
    flirt.run() 
    
    # apply registration to fmap2
    applyxfm = fsl.preprocess.ApplyXFM()
    applyxfm.inputs.in_file = file_fmap_phase
    applyxfm.inputs.reference = os.path.join(path_udata, "median_"+name_udata)
    applyxfm.inputs.in_matrix_file = os.path.join(path_fmap0, "fmap2epi.txt")
    applyxfm.inputs.interp = "trilinear"
    applyxfm.inputs.output_type = "NIFTI"
    applyxfm.inputs.out_file = os.path.join(path_fmap1, "r"+name_fmap1)
    applyxfm.inputs.apply_xfm = True
    applyxfm.run() 
    
    # apply skullstrip mask to fmap1 and fmap2 and save with same header information
    fmap1_img = nb.load(os.path.join(path_fmap0, "r"+name_fmap0))
    arr_fmap1 = fmap1_img.get_fdata()
    fmap2_img = nb.load(os.path.join(path_fmap1, "r"+name_fmap1))
    arr_fmap2 = fmap2_img.get_fdata()
    mask_img = nb.load(os.path.join(path_udata, "mask_median_"+name_udata))
    arr_mask = mask_img.get_fdata()
    
    arr_fmap1 = arr_fmap1 * arr_mask
    arr_fmap2 = (arr_fmap2 * arr_mask) 
    arr_fmap2 = arr_fmap2 + np.abs(np.min(arr_fmap2))
    arr_fmap2 = arr_fmap2 / np.max(arr_fmap2) * 4095 # rescale phase image to be within 0-4095
    
    fmap1_img = nb.Nifti1Image(arr_fmap1, fmap1_img.affine, fmap1_img.header)
    nb.save(fmap1_img, os.path.join(path_fmap0, "pr"+name_fmap0))
    fmap2_img = nb.Nifti1Image(arr_fmap2, fmap1_img.affine, fmap1_img.header)
    nb.save(fmap2_img, os.path.join(path_fmap1, "pr"+name_fmap1))
    
    # prepare fieldmap (saves fieldmap in rad/s)
    prepare = fsl.PrepareFieldmap()
    prepare.inputs.in_magnitude = os.path.join(path_fmap0, "pr"+name_fmap0)
    prepare.inputs.in_phase = os.path.join(path_fmap1, "pr"+name_fmap1)
    prepare.inputs.out_fieldmap = os.path.join(path_fmap0, "fieldmap.nii")
    prepare.inputs.delta_TE = delta_te
    prepare.inputs.scanner = "SIEMENS"
    prepare.inputs.output_type = "NIFTI" 
    prepare.run() 
    
    # effective echo spacing in s
    dwell_time = 1/(bw * ImageMatrixPhaseEncode)
    
    # unmask fieldmap (fsl.FUGUE)
    fugue = fsl.preprocess.FUGUE()
    fugue.inputs.in_file = os.path.join(path_udata, name_udata)
    fugue.inputs.dwell_time = dwell_time
    fugue.inputs.fmap_in_file = os.path.join(path_fmap0, "fieldmap.nii")
    fugue.inputs.smooth3d = smooth
    fugue.inputs.unwarp_direction = udir
    fugue.inputs.save_shift = True
    fugue.inputs.shift_out_file = os.path.join(path_fmap0, "vdm.nii")
    fugue.inputs.output_type = "NIFTI"
    fugue.run() 
    
    # warp coordinate mapping
    generate_coordinate_mapping(file_epi, 
                                0, 
                                path_fmap0, 
                                suffix="fmap", 
                                time=False, 
                                write_output=True)
    
    # apply inverse fieldmap to coordinate mapping
    fugue = fsl.preprocess.FUGUE()
    fugue.inputs.in_file = os.path.join(path_fmap0, "cmap_fmap.nii")
    fugue.inputs.shift_in_file = os.path.join(path_fmap0, "vdm.nii")
    fugue.inputs.forward_warping = False
    fugue.inputs.unwarp_direction = udir
    fugue.inputs.output_type = "NIFTI"
    fugue.run()
    
    # apply cmap to surface
    for i in range(len(file_surf)): 
        path_surf, hemi, name_surf = get_filename(file_surf[i])
        deform_surface(input_surf=file_surf[i], 
                       input_orig=os.path.join(path_udata, "median_"+name_udata), 
                       input_deform=os.path.join(path_fmap0, "cmap_fmap_unwarped.nii"), 
                       input_target=os.path.join(path_udata, "median_"+name_udata), 
                       hemi=hemi, 
                       path_output=path_surf, 
                       interp_method="nearest",
                       smooth_iter=10, 
                       sort_faces=False, 
                       flip_faces=False, 
                       cleanup=True)
    
    # delete created files
    if cleanup:
        os.remove(os.path.join(path_fmap0, "cmap_fmap.nii"))
        os.remove(os.path.join(path_fmap0, "cmap_fmap_unwarped.nii"))
        os.remove(os.path.join(path_fmap0, "fieldmap.nii"))
        os.remove(os.path.join(path_fmap0, "fmap2epi.txt"))
        os.remove(os.path.join(path_fmap1, os.path.splitext(name_fmap1)[0]+"_flirt.mat"))
        os.remove(os.path.join(path_fmap0, "r"+name_fmap0))
        os.remove(os.path.join(path_fmap0, "pr"+name_fmap0))
        os.remove(os.path.join(path_fmap1, "r"+name_fmap1))
        os.remove(os.path.join(path_fmap1, "pr"+name_fmap1))
        os.remove(os.path.join(path_fmap0, os.path.splitext(name_udata)[0])+"_unwarped.nii")
        os.remove(os.path.join(path_fmap0, "vdm.nii"))
        os.remove(os.path.join(path_udata, "mask_median_"+name_udata))
        os.remove(os.path.join(path_udata, "median_"+name_udata))
        os.remove(os.path.join(path_udata, "pmedian_"+name_udata))