def get_weighted_vfs(input_vfs, input_snr, hemi, path_output):
    """
    This function scales the visual fieldsign assignments using a given snr estimate.
    Inputs:
        *input_vfs: input file with visual fieldsign assignments.
        *input_snr: input file with snr estimates.
        *hemi: hemisphere.
        *path_output: output folder where output is saved.
    
    created by Daniel Haenelt
    Date created: 13-02-2019
    Last modified: 13-02-2019
    """
    import os
    import numpy as np
    import nibabel as nb

    # load visual fieldsign map
    vfs_img = nb.load(input_vfs)
    vfs_array = vfs_img.get_fdata()

    # load snr map
    snr_img = nb.load(input_snr)
    snr_array = snr_img.get_fdata()
    
    # normalize snr map
    snr_array = snr_array / np.max(snr_array)    

    # get weighted visual fieldsign map
    vfs_array = vfs_array * snr_array

    # write output
    out_img = nb.Nifti1Image(vfs_array, vfs_img.affine, vfs_img.header)
    nb.save(out_img, os.path.join(path_output,hemi+".fieldsign_weighted.mgh"))