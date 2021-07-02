# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import nibabel as nb


def include_pial_correction(path, sub):
    """Include pial correction.

    This function takes manual corrections from the file pial_edit.mgz which 
    shall be located in the freesurfer mri folder and includes them in the 
    brainmask. The manual corrected brainmask is saved as 
    brain.finalsurfs.manedit.mgz. Note that the corrections in the pial_edit.mgz 
    are done with brush value 256 and eraser value -1.    

    Parameters
    ----------
    path : str
        Path to SUBJECTS_ID.
    sub : str
        Freesurfer subject name.

    Returns
    -------
    None.
    
    """

    # open pial edit
    edit_img = nb.load(os.path.join(path, sub, "mri", "pial_edit.mgz"))
    edit_array = edit_img.get_fdata()

    # open brainmask
    brainmask_img = nb.load(os.path.join(path, sub, "mri", "brainmask.mgz"))
    brainmask_array = brainmask_img.get_fdata()

    # include manual edits to brainmask
    brainmask_array[edit_array == 256] = 255
    brainmask_array[edit_array == -1] = 1

    # save brainmask as brain.finalsurfs.manedit.mgz
    output = nb.Nifti1Image(brainmask_array, brainmask_img.affine,
                            brainmask_img.header)
    nb.save(output,
            os.path.join(path, sub, "mri", "brain.finalsurfs.manedit.mgz"))
