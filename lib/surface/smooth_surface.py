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
    import sys
    import subprocess
    from lib.io.get_filename import get_filename
       
    # make output folder
    path_output, _, _ = get_filename(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # smooth surface
    try:
        subprocess.run(['mris_smooth', '-n', str(n_iter), '-nw', file_in, file_out], check = True)
    except subprocess.CalledProcessError:
        sys.exit("Surface smoothing failed!")