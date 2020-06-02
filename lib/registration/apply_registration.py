def apply_registration(file_in, cmap_in, file_out, interpolation="linear", r=[0.4,0.4,0.4]):
    """
    This function applies a coordinate mapping to a volume. Optionally, the voxel size of the output
    volume can be changed. This is achieved by adjusting the coordinate mapping to the new voxel 
    size before application.
    Inputs:
        *file_in: filename of input volume.
        *cmap_in: filename of coordinate mapping.
        *file_out: filename of output volume.
        *interpolation: interpolation type (linear or nearest).
        *r: destination voxel size after upsampling (performed if not None).
    Outputs:
        *nibabel object instance of transformed input.
    
    created by Daniel Haenelt
    Date created: 30-05-2020
    Last modified: 02-06-2020
    """
    import os
    import numpy as np
    import nibabel as nb
    from nighres.registration import apply_coordinate_mappings
    from lib.io.get_filename import get_filename
    from lib.utils.upsample_volume import upsample_volume
       
    # make output folder     
    path_output = os.path.dirname(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)
      
    # filename for temporary cmap copy
    _, _, ext_cmap = get_filename(cmap_in)  
    tmp = np.random.randint(0, 10, 5)
    tmp_string = ''.join(str(i) for i in tmp)
    file_tmp = os.path.join(path_output,"tmp_"+tmp_string+ext_cmap)
    file_tmp2 = os.path.join(path_output,"tmp2_"+tmp_string+ext_cmap)
    
    # adjust coordinate mapping
    if r:
    
        # upsample cmap
        upsample_volume(cmap_in, file_tmp, dxyz = r, rmode = "Linear")
        upsample_volume(cmap_in, file_tmp2, dxyz = r, rmode = "NN")
        
        # mask upsampled cmap
        cmap = nb.load(file_tmp)
        mask = nb.load(file_tmp2)
        
        cmap_array = cmap.get_fdata()
        mask_array = mask.get_fdata()
            
        mask_array = np.sum(mask_array, axis=3)
        mask_array[mask_array != 0] = 1
        
        cmap_array[:,:,:,0][mask_array == 0] = 0
        cmap_array[:,:,:,1][mask_array == 0] = 0
        cmap_array[:,:,:,2][mask_array == 0] = 0
            
        cmap = nb.Nifti1Image(cmap_array, cmap.affine, cmap.header)
    
    else:
        
        cmap = nb.load(cmap_in)
    
    # apply coordinate mapping
    res = apply_coordinate_mappings(image = file_in, # input 
                                    mapping1 = cmap, # cmap
                                    interpolation = interpolation, # nearest or linear
                                    padding = "zero", # closest, zero or max
                                    save_data = False, # save output data to file (boolean)
                                    overwrite = False, # overwrite existing results (boolean)
                                    output_dir = None, # output directory
                                    file_name = None, # base name with file extension for output
                                    )
    
    # write output
    nb.save(res["result"], file_out)
    
    # remove temporary files
    if r:
        os.remove(file_tmp)
        os.remove(file_tmp2)

    return res["result"]