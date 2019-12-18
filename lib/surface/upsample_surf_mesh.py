def upsample_surf_mesh(file_in, file_out, niter, method):
    """
    The scripts takes generated FreeSurfer surfaces and upsamples them using Jon Polimeni's function 
    mris_mesh_subdivide.
    Inputs:
        *file_in: filename of input surface.
        *file_out: filename of output surface.
        *niter: number of upsampling iterations.
        *method: upsampling method (linear, loop, butterfly).

    created by Daniel Haenelt
    Date created: 01-11-2018             
    Last modified: 18-12-2019
    """
    import os

    # make output folder
    path_output = os.path.dirname(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # subdivide the surfaces
    os.system("mris_mesh_subdivide" + \
              " --surf " + file_in + \
              " --out " + file_out + \
              " --method " + method + \
              " --iter " + str(niter))