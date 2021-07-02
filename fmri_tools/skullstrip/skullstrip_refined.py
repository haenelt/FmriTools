# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import nibabel as nb
from nipype.interfaces.freesurfer import ApplyVolTransform

# local inputs
from ..io.get_filename import get_filename
    
    
def skullstrip_refined(file_mask1, file_mask2):
    """Skullstrip refined.

    The purpose of the following function is to enhance the skullstrip mask in 
    native space. It uses a second mask which was manually corrected during the 
    freesurfer segmentation. This corrected brainmask is converted to native 
    space and multiplied with the initial brainmask.    

    Parameters
    ----------
    file_mask1 : str
        Brainmask in original space.
    file_mask2 : str
        Manually corrected brainmask in freesurfer space (brain.finalsurfs.mgz).

    Returns
    -------
    file_out : str
        Filename of enhanced brainmask.
    
    """
    
    # get output path and basename
    path_output, name_output, _ = get_filename(file_mask1)
    
    # filename of temporary and enhanced brainmask
    file_temp = os.path.join(path_output, "temp.nii")
    file_out = os.path.join(path_output, name_output+"_enhanced.nii")

    # bring skullstrip_mask from conformed space into original space
    transmask = ApplyVolTransform()
    transmask.inputs.source_file = file_mask2
    transmask.inputs.target_file = file_mask1
    transmask.inputs.reg_header = True
    transmask.inputs.interp = "nearest"
    transmask.inputs.transformed_file = file_temp
    transmask.inputs.args = "--no-save-reg"
    transmask.run()

    # load first brainmask in original space
    mask1 = nb.load(file_mask1)
    mask1_array = mask1.get_fdata()

    # load second brainmask transformed into original space
    mask2 = nb.load(file_temp)
    mask2_array = mask2.get_fdata()

    # make second brainmask binary
    mask2_array[mask2_array == 1] = 0
    mask2_array[mask2_array > 0] = 1

    # multiply both masks
    mask_enhanced_array = mask1_array * mask2_array
    
    # write enhancec brainmask
    mask_enhanced = nb.Nifti1Image(mask_enhanced_array, mask1.affine,
                                   mask1.header)
    nb.save(mask_enhanced, file_out)

    # remove temporary file
    os.remove(file_temp)
    
    return file_out
