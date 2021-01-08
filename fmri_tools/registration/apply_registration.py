# -*- coding: utf-8 -*-

# python standard library inputs
import os
import datetime

# external inputs
import numpy as np
import nibabel as nb
from nighres.registration import apply_coordinate_mappings

# local inputs
from fmri_tools.io.get_filename import get_filename
from fmri_tools.utils.resample_volume import resample_volume


def apply_registration(file_in, cmap_in, file_out, interpolation="linear", 
                       r=[0.4,0.4,0.4]):
    """ Apply registration

    This function applies a coordinate mapping to a volume. Optionally, the 
    voxel size of the output volume can be changed. This is achieved by 
    adjusting the coordinate mapping to the new voxel size before application.    

    Parameters
    ----------
    file_in : str
        Filename of input volume.
    cmap_in : str
        Filename of coordinate mapping.
    file_out : str
        Filename of output volume.
    interpolation : str, optional
        Interpolation type (linear or nearest). The default is "linear".
    r : list, optional
        Destination voxel size after upsampling (performed if not None). The 
        default is [0.4,0.4,0.4].

    Returns
    -------
    niimg
        Transformed volume.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 30-05-2020
    Last modified: 25-10-2020
    
    """
          
    # make output folder     
    path_output = os.path.dirname(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)
      
    # filename for temporary cmap copy
    _, _, ext_cmap = get_filename(cmap_in)      
    tmp1 = np.random.randint(0, 10, 5)
    tmp1 = ''.join(str(i) for i in tmp1)
    tmp2 = datetime.datetime.now().strftime("%S%f")
    tmp_string = tmp1 + tmp2
    file_tmp = os.path.join(path_output,"tmp_"+tmp_string+ext_cmap)
    file_tmp2 = os.path.join(path_output,"tmp2_"+tmp_string+ext_cmap)
    
    if os.path.exists(file_tmp) or os.path.exists(file_tmp2):
        raise FileExistsError("Temporary file already exists!")
    
    # adjust coordinate mapping
    if r:
    
        # resample cmap
        resample_volume(cmap_in, file_tmp, dxyz = r, rmode = "Linear")
        resample_volume(cmap_in, file_tmp2, dxyz = r, rmode = "NN")
        
        # mask resampled cmap
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
