# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
from nipype.interfaces.freesurfer.preprocess import MRIConvert

# local inputs
from ..io.get_filename import get_filename


def mgh2nii(file_in, path_output, out_type="nii"):
    """MGH2NII.

    This function converts a volume file from freesurfer mgh to nifti format.    

    Parameters
    ----------
    file_in : str
        Full path of the input file.
    path_output : str
        Path where output is written.
    out_type : str, optional
        Target type of file. The default is "nii".

    Returns
    -------
    None.
    
    """

    # get filename
    path, name, ext = get_filename(file_in)

    # convert volume to nifti format
    mc = MRIConvert()
    mc.inputs.in_file = file_in
    mc.inputs.out_file = os.path.join(path_output, name + "." + out_type)
    mc.inputs.in_type = ext.replace('.', '')
    mc.inputs.out_type = out_type
    mc.run()
