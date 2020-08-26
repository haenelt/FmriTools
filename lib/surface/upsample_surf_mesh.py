def upsample_surf_mesh(file_in, file_out, n_iter, method):
    """
    The scripts takes generated FreeSurfer surfaces and upsamples them using Jon Polimeni's function 
    mris_mesh_subdivide.
    Inputs:
        *file_in: filename of input surface.
        *file_out: filename of output surface.
        *n_iter: number of upsampling iterations.
        *method: upsampling method (linear, loop, butterfly).

    created by Daniel Haenelt
    Date created: 01-11-2018             
    Last modified: 26-08-2020
    """
    import os
    import sys
    import subprocess
    from lib.io.get_filename import get_filename

    # make output folder
    path_output, _, _ = get_filename(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # subdivide surface
    try:
        subprocess.run(['mris_mesh_subdivide', 
                        '--surf', file_in, 
                        '--out', file_out, 
                        '--method', method, 
                        '--iter', str(n_iter)], check = True)
    except subprocess.CalledProcessError:
        sys.exit("Surface subdivision failed!")