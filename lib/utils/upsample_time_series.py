def upsample_time_series(input, n=2):
    """
    This function upsamples a time series using the afni function 3dUpsample. Before running the
    function, set the afni environment by calling AFNI in the terminal. Output is an upsampled nifti
    time series.
    Inputs:
        *input: time series filename.
        *n: upsampling factor.
        
    created by Daniel Haenelt
    Date created: 20-09-2019           
    Last modified: 20-09-2019
    """
    import os
    from sh import gunzip
    
    clean_unzip = 0
    if os.path.splitext(input)[1] == ".gz":
        gunzip(input)
        clean_unzip = 1
        input = os.path.splitext(input)[0]
        
    # prepare path and filename
    path_file = os.path.dirname(input)
    name_file = os.path.splotext(os.path.basename(input))[0]

    # upsample vaso and bold time series
    os.system("3dUpsample -overwrite -datum short " + \
              "-prefix " + os.path.join(path_file,name_file + "_upsampled.nii") + \
              " -n " + str(n) + " -input " + input)
    
    if clean_unzip:
        os.remove(input)