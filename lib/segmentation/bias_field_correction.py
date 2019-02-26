def bias_field_correction(filename, pathSPM12):
    """
    Renzo recommended to perform a bias field correction using SPM12 before doing the segmentation 
    with FreeSurfer. FreeSurfer recommends Bias FWHM = 18 and Sampling distance = 2 for MEMPRAGE at 
    7 T, which is also set here. Outputs are saved in the input folder.
    Inputs:
        *filename: path of input image.
        *prefix: Defined prefix for the output image.
        *pathSPM12: path to SPM12 toolbox.
        
    created by Daniel Haenelt
    Date created: 01-11-2018             
    Last modified: 17-12-2018
    """
    import os
    from nipype.interfaces.spm import NewSegment
    from nipype.interfaces.matlab import MatlabCommand
    
    # set matlab path to SPM12 folder
    MatlabCommand.set_default_paths(pathSPM12)

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