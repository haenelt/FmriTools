def make_sphere(file_in, file_out, n_inflate=100, radius=None):
    """
    The scripts takes a generated FreeSurfer mesh and transformes it into
    a sphere with defined radius.
    Inputs:
        *file_in: filename of input surface.
        *file_out: filename of output surface.
        *n_inflated: number of inflating iterations (if > 0).
        *radius: radius of final sphere in mm (if not None).

    created by Daniel Haenelt
    Date created: 26-08-2020       
    Last modified: 26-08-2020
    """
    import os
    import sys
    import subprocess
    import numpy as np
    from shutil import copyfile
    from nibabel.freesurfer.io import read_geometry, write_geometry
    from lib.io.get_filename import get_filename
    from lib.surface.inflate_surf_mesh import inflate_surf_mesh
    
    def cart2pol(x, y, z):
        r = np.sqrt(x**2+y**2+z**2)
        phi = np.arctan2(y, x)
        theta = np.arccos(z/r)
        return r, phi, theta
    
    def pol2cart(r, phi, theta):
        x = r*np.sin(theta)*np.cos(phi)
        y = r*np.sin(theta)*np.sin(phi)
        z = r*np.cos(theta)
        return x, y, z

    # make output folder
    path_output, _, _ = get_filename(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)
    
    # temporary file
    tmp = np.random.randint(0, 10, 5)
    tmp_string = ''.join(str(i) for i in tmp)
    file_tmp = os.path.join(path_output, tmp_string)

    # inflate surface mesh
    if n_inflate:
        inflate_surf_mesh(file_in, 
                          file_tmp, 
                          n_inflate)
    else:
        copyfile(file_in, file_tmp)
        
    
    # inflate surface
    try:
        subprocess.run(['mris_sphere', 
                        '-q', 
                        file_tmp, 
                        file_out], check = True)
    except subprocess.CalledProcessError:
        sys.exit("Sphere computation failed!")
    

    # change radius      
    if radius:
        vtx, fac = read_geometry(file_out)
        r, phi, theta = cart2pol(vtx[:,0], vtx[:,1], vtx[:,2])
        r[:] = radius
        vtx[:,0], vtx[:,1], vtx[:,2] = pol2cart(r, phi, theta)
        write_geometry(file_out, vtx, fac)
    
    # remove temporary file
    os.remove(file_tmp)
    