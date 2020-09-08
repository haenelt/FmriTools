def write_label(arr_label, file_out):
    """
    This function writes a textfile which can be read as label file in 
    freesurfer.
    Inputs:
        *arr_label: list with label indices.
        *file_out: filename of label file.

    created by Daniel Haenelt
    Date created: 03-09-2020
    Last modified: 03-09-2020
    """
    
    # number of labels
    n_label = len(arr_label)
    
    with open(file_out, "w") as f:
        f.write('#!ascii label  , from subject  vox2ras=TkReg\n')
        f.write(str(n_label)+'\n')
        for i in range(n_label):
            f.write(str(arr_label[i])+' 0.000 0.000 0.000 0.000\n')
