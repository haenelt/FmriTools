# -*- coding: utf-8 -*-


def write_label(arr_label, file_out):
    """ Write label

    This function writes a textfile which can be read as label file in 
    freesurfer.        

    Parameters
    ----------
    arr_label : list
        List of label indices.
    file_out : str
        Filename of label file.

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
    Last modified: 20-10-2020

    """
    
    # check filename
    if isinstance(file_out, str):
        if not file_out.endswith("txt"):            
            raise ValueError("Currently supported file formats is txt.")
    else:
        raise ValueError("Filename must be a string!")
    
    # number of labels
    n_label = len(arr_label)
    
    with open(file_out, "w") as f:
        f.write('#!ascii label  , from subject  vox2ras=TkReg\n')
        f.write(str(n_label)+'\n')
        for i in range(n_label):
            f.write(str(arr_label[i])+' 0.000 0.000 0.000 0.000\n')
