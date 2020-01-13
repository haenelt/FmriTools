def upsample_volume(file_in, file_out, dxyz=[0.4, 0.4, 0.4], rmode="Cu"):
    """
    This function upsamples a nifti volume using the afni function 3dresample. Before running the
    function, set the afni environment by calling AFNI in the terminal. Output is an upsampled nifti
    volume.
    Inputs:
        *file_in: nifti input filename.
        *file_out: nifti output filename.
        *dxyz: array of target resolution in single dimensions.
        *rmode: interpolation methods (Linear, NN, Cu, Bk).
        
    created by Daniel Haenelt
    Date created: 16-12-2019        
    Last modified: 16-12-2019
    """
    import os
    from sh import gunzip
    
    clean_unzip = 0
    if os.path.splitext(file_in)[1] == ".gz":
        gunzip(file_in)
        clean_unzip = 1
        file_in = os.path.splitext(file_in)[0]

    # upsample volume
    os.system("3dresample " + \
              "-dxyz " + str(dxyz[0]) + " " + str(dxyz[1]) + " " + str(dxyz[2]) + " " +\
              "-rmode " + str(rmode) + " " + \
              "-inset " + file_in + " " + \
              "-prefix " + file_out)
    
    if clean_unzip:
        os.remove(file_in)