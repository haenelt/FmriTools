def clean_coordinate_mapping(cmap_source, cmap_target, overwrite_file=True, save_mask=False):
    """
    Voxels in the target coordinate mapping are masked out based on found voxel displacements in the
    source coordinate mapping. This is done to remove smeared regions caused by interpolations with
    background values in the case of deforming a slab within a larger image array.
    Inputs:
        *cmap_source: filename of source coordinate mapping.
        *cmap_target: filename of target coordinate mapping.
        *overwrite_file: overwrite target coordinate mapping (boolean).
        *save_mask: write out mask (boolean).
    Outputs:
        *results: nibabel instances of cleaned cmap and corresponding mask (dict).

    created by Daniel Haenelt
    Date created: 23-05-2020             
    Last modified: 23-05-2020
    """
    import os
    import numpy as np
    import nibabel as nb
    from lib.io.get_filename import get_filename
    
    # get filename    
    path_file, _, _ = get_filename(cmap_target)
    
    # load data
    cmap1_img = nb.load(cmap_source)
    cmap1_array = cmap1_img.get_fdata()
    
    cmap2_img = nb.load(cmap_target)
    cmap2_array = cmap2_img.get_fdata()
    
    mask_img = nb.load(cmap_target)
    mask_img.header["dim"][0] = 3
    mask_img.header["dim"][4] = 1
    mask_array = np.zeros_like(mask_img.get_fdata()[:,:,:,0])
    
    x_max = cmap2_img.header["dim"][1]
    y_max = cmap2_img.header["dim"][2]
    z_max = cmap2_img.header["dim"][3]
    
    # get nearest voxel coordinates
    x0 = np.floor(cmap1_array[:,:,:,0].flatten()).astype(int)
    x1 = np.ceil(cmap1_array[:,:,:,0].flatten()).astype(int)
    y0 = np.floor(cmap1_array[:,:,:,1].flatten()).astype(int)
    y1 = np.ceil(cmap1_array[:,:,:,1].flatten()).astype(int)
    z0 = np.floor(cmap1_array[:,:,:,2].flatten()).astype(int)
    z1 = np.ceil(cmap1_array[:,:,:,2].flatten()).astype(int)
    
    # exclude voxels which do not fit in the target array
    outlier = []
    outlier.extend(np.where(x0 < 0)[0])
    outlier.extend(np.where(x1 < 0)[0])
    outlier.extend(np.where(y0 < 0)[0])
    outlier.extend(np.where(y1 < 0)[0])
    outlier.extend(np.where(z0 < 0)[0])
    outlier.extend(np.where(z1 < 0)[0])
    outlier.extend(np.where(x0 >= x_max)[0])
    outlier.extend(np.where(x1 >= x_max)[0])
    outlier.extend(np.where(y0 >= y_max)[0])
    outlier.extend(np.where(y1 >= y_max)[0])
    outlier.extend(np.where(z0 >= z_max)[0])
    outlier.extend(np.where(z1 >= z_max)[0])
    outlier = list(np.unique(outlier))
    
    x0 = np.delete(x0, outlier)
    x1 = np.delete(x1, outlier)
    y0 = np.delete(y0, outlier)
    y1 = np.delete(y1, outlier)
    z0 = np.delete(z0, outlier)
    z1 = np.delete(z1, outlier)
    
    # get final mask
    mask_array[x0,y0,z0] = 1
    mask_array[x1,y1,z1] = 1
    mask_array[x1,y0,z0] = 1
    mask_array[x0,y1,z0] = 1
    mask_array[x0,y0,z1] = 1
    mask_array[x1,y1,z0] = 1
    mask_array[x0,y1,z1] = 1
    mask_array[x1,y0,z1] = 1
    
    # apply mask to cmap
    cmap2_array[:,:,:,0] *= mask_array
    cmap2_array[:,:,:,1] *= mask_array
    cmap2_array[:,:,:,2] *= mask_array
       
    # get output
    results = dict()
    results["cmap"] = nb.Nifti1Image(cmap2_array, cmap2_img.affine, cmap2_img.header)
    results["mask"] = nb.Nifti1Image(mask_array, mask_img.affine, mask_img.header)
    
    # write output
    if overwrite_file:
        nb.save(results["cmap"], cmap_target)
    
    if save_mask:
        nb.save(results["mask"], os.path.join(path_file,"cmap_mask.nii"))
    
    return results
