def inflate_surf_mesh(file_in, file_out, n_iter):
    """
    The scripts takes a generated FreeSurfer surfaces and inflates it.
    Inputs:
        *file_in: filename of input surface.
        *file_out: filename of output surface.
        *n_iter: number of inflating iterations.

    created by Daniel Haenelt
    Date created: 17-12-2019             
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

    # inflate surface
    try:
        subprocess.run(['mris_inflate', 
                        '-n', str(n_iter), 
                        '-no-save-sulc', 
                        file_in, 
                        file_out], check = True)
    except subprocess.CalledProcessError:
        sys.exit("Surface inflation failed!")