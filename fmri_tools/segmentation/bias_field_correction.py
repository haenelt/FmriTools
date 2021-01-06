# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
from nipype.interfaces.spm import NewSegment


def bias_field_correction(filename):
    """ Bias field correction

    Renzo recommended to perform a bias field correction using SPM12 before 
    doing the segmentation with FreeSurfer. FreeSurfer recommends Bias FWHM = 18 
    and Sampling distance = 2 for MEMPRAGE at 7 T, which is also set here. 
    Outputs are saved in the input folder.    

    Parameters
    ----------
    filename : str
        Path of input image.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 01-11-2018             
    Last modified: 12-10-2020
    
    """
       
    # get path of filename
    path = os.path.dirname(filename)

    # bias field correction
    os.chdir(path)
    bias = NewSegment()
    bias.inputs.channel_files = os.path.join(filename)
    bias.inputs.channel_info = (0.001, 18, (True, True))
    bias.inputs.affine_regularization = "mni"
    bias.inputs.sampling_distance = 2
    bias.inputs.use_v8struct = True
    bias.inputs.warping_regularization = [0, 0.001, 0.5, 0.05, 0.2]
    bias.inputs.write_deformation_fields = [False, False]
    bias.inputs.mfile = True
    tissue1 = ((os.path.join(pathSPM12, "tpm/TPM.nii"), 1), 1, (False,False), (False, False))
    tissue2 = ((os.path.join(pathSPM12, "tpm/TPM.nii"), 2), 1, (False,False), (False, False))
    tissue3 = ((os.path.join(pathSPM12, "tpm/TPM.nii"), 3), 2, (False,False), (False, False))
    tissue4 = ((os.path.join(pathSPM12, "tpm/TPM.nii"), 4), 3, (False,False), (False, False))
    tissue5 = ((os.path.join(pathSPM12, "tpm/TPM.nii"), 5), 4, (False,False), (False, False))
    tissue6 = ((os.path.join(pathSPM12, "tpm/TPM.nii"), 6), 2, (False,False), (False, False))
    bias.inputs.tissues = [tissue1, tissue2, tissue3, tissue4, tissue5, tissue6]
    bias.run()
    