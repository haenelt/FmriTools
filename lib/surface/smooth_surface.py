def smooth_surface(file_in, file_out, n_iter):
    """
    This function smoothes a surface file using freesurfer.
    Inputs:
        *file_in: filename of input surface.
        *file_out: filename of output surface.
        *n_iter: number of smoothing iterations.
        
    created by Daniel Haenelt
    Date created: 13-07-2019
    Last modified: 18-12-2019
    """
    from nipype.interfaces.freesurfer import SmoothTessellation
         
    # smooth surface
    smooth = SmoothTessellation()
    smooth.inputs.in_file = file_in
    smooth.inputs.out_file = file_out
    smooth.inputs.smoothing_iterations = n_iter
    smooth.inputs.disable_estimates = True
    smooth.run()