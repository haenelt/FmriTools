def average_layer(img_input, path_output, basename_output):
    """
    This averages data across different layers. Input arrays should be in mgh format. The output 
    gets the suffix _avg.
    Inputs:
        *img_input: list of input layers.
        *path_output: path where output is saved.
        *basename_output: basename of written output file.
    
    created by Daniel Haenelt
    Date created: 25-10-2019
    Last modified: 25-10-2019
    """
    import os
    import numpy as np
    import nibabel as nb
    
    # initialise array
    data = nb.load(img_input[0])
    data_array = np.zeros_like(data.get_fdata())
    
    # sum over layers
    for i in range(len(img_input)):
        data_array += nb.load(img_input[i]).get_fdata()    

    # divide by number of layers
    data_array = data_array / len(img_input)
    
    # write output file
    output = nb.Nifti1Image(data_array,data.affine,data.header)
    nb.save(output,os.path.join(path_output,basename_output+"_avg.mgh"))