def morph2dense(source_sphere,target_sphere,input_morph,path_output):
    """
    This function maps a morphological file from a source surface to a target target surface.    
    Inputs:
        *source_sphere: source surface.
        *target_sphere: target surface.
        *input_morph: morphological input file.
        *path_output: path where output is saved.
        
    created by Daniel Haenelt
    Date created: 13-07-2019
    Last modified: 13-07-2019
    """
    import os
    from nibabel.freesurfer.io import read_morph_data, write_morph_data, read_geometry
    from scipy.interpolate import griddata    

    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # transform morphological data to dense surfaces
    pts_sphere_dense, _ = read_geometry(target_sphere)
    pts_sphere, _ = read_geometry(source_sphere)
    
    # get morphological data
    morph = read_morph_data(input_morph)
    
    # do the transformation
    method = "nearest"
    morph_dense = griddata(pts_sphere, morph, pts_sphere_dense, method)
        
    # write dense morphological data
    write_morph_data(os.path.join(path_output,os.path.basename(input_morph)), morph_dense)