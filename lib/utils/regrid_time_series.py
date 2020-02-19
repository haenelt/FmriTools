def regrid_time_series(input, path_output, TR_source, TR_target, t_start=0):
    """
    This function interpolates the time series onto a new time grid using cubic interpolation.
    Inputs:
        *input: time series filename.
        *n: upsampling factor.
        
    created by Daniel Haenelt
    Date created: 19-02-2020           
    Last modified: 19-02-2020
    """
    import os
    import numpy as np
    import nibabel as nb
    from scipy.interpolate import griddata
    from lib.io.get_filename import get_filename

    # get filename
    _, name_input, ext_input = get_filename(input)    
    
    # load input time series
    data = nb.load(input)
    data_array = data.get_fdata()

    # get matrix dimensions
    nx = np.shape(data_array)[0]
    ny = np.shape(data_array)[1]
    nz = np.shape(data_array)[2]
    nt = np.shape(data_array)[3]

    # source time points
    t = TR_source * np.arange(0,nt) + t_start

    # target time points
    t_end = np.max(t) + 10
    t_new = np.arange(0,t_end, TR_target)
    t_new = t_new[t_new <= np.max(t)]

    data_array_new = np.zeros((nx,ny,nz,len(t_new)))
    for x in range(nx):    
        for y in range(ny):
            for z in range(nz):
                data_array_new[x,y,z,:] = griddata(t, data_array[x,y,z,:], t_new, method="cubic")
    
    # clean upsampled array
    data_array_new[np.isnan(data_array_new)] = 0
    data_array_new[data_array_new < 0] = 0
    
    # update data header
    data.header["dim"][4] = len(t_new)
    data.header["datatype"] = 16

    # write output
    output = nb.Nifti1Image(data_array_new, data.affine, data.header)
    nb.save(output, os.path.join(path_output,name_input+"_upsampled"+ext_input))
    
    # change TR in header
    os.system("3drefit " + \
              "-TR " + str(TR_target) + " " + \
              os.path.join(path_output,name_input+"_upsampled"+ext_input))

def regrid_time_series_afni(input, n=2):
    """
    This function upsamples a time series using the afni function 3dUpsample. Before running the
    function, set the afni environment by calling AFNI in the terminal. Output is an upsampled nifti
    time series.
    Inputs:
        *input: time series filename.
        *n: upsampling factor.
        
    created by Daniel Haenelt
    Date created: 20-09-2019           
    Last modified: 19-02-2020
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
    name_file = os.path.splitext(os.path.basename(input))[0]

    # upsample vaso and bold time series
    os.system("3dUpsample -overwrite -datum short " + \
              "-prefix " + os.path.join(path_file,name_file + "_upsampled.nii") + \
              " -n " + str(n) + " -input " + input)
    
    if clean_unzip:
        os.remove(input)