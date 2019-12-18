def get_curvature(file_in, path_output, a=10):
    """
    This function calculates a curvature file for an input surface mesh using freesurfer. The input 
    file needs to have a prefix which indicates the hemisphere of the surface mesh.
    Inputs:
        *path: filename of input surface.
        *path_output: path where output is written.
        *a: number of smoothing iterations.
        
    created by Daniel Haenelt
    Date created: 13-07-2019
    Last modified: 18-12-2019
    """
    import sys
    import os
    from nipype.interfaces.freesurfer import Curvature
    
    # get hemi from filename
    hemi = os.path.splitext(os.path.basename(file_in))[0]
    if not hemi == "lh" or hemi == "rh":
        sys.exit("Could not identify hemi from filename!")
    
    # calculate curvature file
    curv = Curvature()
    curv.inputs.in_file = file_in
    curv.inputs.save = True
    curv.inputs.averages = a # for curv file
    curv.run()
    
    # rename mean curvature to curv
    os.rename(file_in+".H", os.path.join(path_output,hemi+".curv"))
    
    # delete gaussian curvature file and temporary input file
    os.remove(file_in+".K")