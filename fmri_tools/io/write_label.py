# -*- coding: utf-8 -*-

# python standard library inputs
import os

# local inputs
from fmri_tools.io.get_filename import get_filename


def write_label(file_out, arr_label):
    """ Write label

    This function writes a textfile which can be read as label file in 
    freesurfer.        

    Parameters
    ----------
    file_out : str
        Filename of label file.
    arr_label : list
        List of label indices.

    Raises
    ------
    ValueError
        If `file_out` is not a string or has a file extension which is not 
        supported.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 03-09-2020
    Last modified: 19-11-2020

    """
    
    # check filename
    if isinstance(file_out, str):
        if not file_out.endswith("label"):            
            raise ValueError("Currently supported file format is txt.")
    else:
        raise ValueError("Filename must be a string!")
    
    # make output folder
    path_output, _, _ = get_filename(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)
    
    # number of labels
    n_label = len(arr_label)
    
    with open(file_out, "w") as f:
        f.write('#!ascii label  , from subject  vox2ras=TkReg\n')
        f.write(str(n_label)+'\n')
        for i in range(n_label):
            f.write(str(arr_label[i])+' 0.000 0.000 0.000 0.000\n')
