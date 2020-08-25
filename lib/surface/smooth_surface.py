def smooth_surface(file_in, file_out, n_iter):
    """
    This function smoothes a surface mesh using freesurfer.
    Inputs:
        *file_in: filename of input surface.
        *file_out: filename of output surface.
        *n_iter: number of smoothing iterations.
        
    created by Daniel Haenelt
    Date created: 13-07-2019
    Last modified: 25-08-2020
    """
    import os
         
    # smooth surface   
    os.system("mris_smooth" + \
              " -n "+str(n_iter) + \
              " " + file_in + \
              " " + file_out)