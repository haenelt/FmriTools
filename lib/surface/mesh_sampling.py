def mesh_sampling(surf_in, file_in, boundaries_in, path_output, layer, r=[0.4,0.4,0.4]):
    """
    This function samples data from an image volume to a surface mesh from specific layers defined 
    by a levelset image.
    Inputs:
        *surf_in: filename of input surface mesh.
        *file_in: filename of input volume from which data is sampled.
        *boundaries_in: filename of 4D levelset image.
        *path_output: path where output is written.
        *layer: which layers to sample (array of integers).
        *r: destination voxel size after upsampling (performed if not None).
    
    created by Daniel Haenelt
    Date created: 18-12-2019
    Last modified: 18-12-2019
    """
    import sys
    import os
    import shutil as sh
    import nibabel as nb
    from os.path import join, exists, basename, splitext
    from nighres.laminar import profile_sampling
    from lib.io.get_filename import get_filename
    from lib.utils.upsample_volume import upsample_volume
    from lib.mapping import map2surface

    # make output folder
    if not exists(path_output):
        os.makedirs(path_output)
    
    # filenames
    _, name_file, ext_file = get_filename(file_in)   
    _, hemi, name_surf = get_filename(surf_in)  
     
    name_surf = name_surf[1:]
    name_profile = splitext(basename(file_in))[0]+"_profile"
    
    # check hemi
    if not hemi == "lh" or hemi == "rh":
        sys.exit("Could not identify hemi from filename!")
    
    # upsample volume
    if not r == None:
        name_file = name_file+"_upsampled"
        upsample_volume(file_in, join(path_output, name_file+ext_file), r, "Cu")
    else:
        if file_in != join(path_output, name_file+ext_file):
            sh.copyfile(file_in, join(path_output, name_file+ext_file))
        
    
    # get profile sampling
    profile = profile_sampling(boundaries_in, 
                               join(path_output, name_file+ext_file),
                               save_data=True, 
                               overwrite=True,
                               output_dir=path_output,
                               file_name="profile")
    
    # rename profile sampling output
    os.rename(join(path_output, "profile_lps-data.nii.gz"),
              join(path_output, name_profile+".nii.gz"))
    
    # map single layers
    for i in range(len(layer)):
        
        # get level
        data = profile["result"]
        data.header["dim"][0] = 1
        data_array = data.get_fdata()
        data_array = data_array[:,:,:,layer[i]]
        out = nb.Nifti1Image(data_array, data.affine, data.header)
        nb.save(out, join(path_output,"temp.nii"))
        
        # do the mapping
        map2surface(surf_in, 
                    join(path_output,"temp.nii"),
                    hemi, 
                    path_output,
                    input_white=None, 
                    input_ind=None, 
                    cleanup=True)
    
        # rename mapping file
        os.rename(join(path_output,hemi+".temp_"+name_surf+"_def.mgh"),
                  join(path_output,hemi+"."+name_file+"_layer"+str(layer[i])+".mgh"))
    
    # clean temp
    os.remove(join(path_output,"temp.nii"))