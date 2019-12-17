def inflate_surf_mesh(file_in, file_out, niter):
    """
    The scripts takes a generated FreeSurfer surfaces and inflates it.
    Inputs:
        *file_in: filename of input surface.
        *file_out: filename of output surface.
        *niter: number of inflating iterations.

    created by Daniel Haenelt
    Date created: 17-12-2019             
    Last modified: 17-12-2019
    """
    import os

    # make output folder
    path_output = os.path.dirname(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # subdivide the surfaces
    os.system("mris_inflate" + \
              " -n " + str(niter) + \
              " -no-save-sulc" + \
              " " + file_in + \
              " " + file_out)
