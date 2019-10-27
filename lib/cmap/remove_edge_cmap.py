def remove_edge_cmap(input_cmap, path_output, basename_output, edge_threshold=2):
    """
    This function removes smeared edges from a coordinate mapping. Depending on the interpolation
    method, coordinate mapping slabs within a larger volume can have blurred edges where voxels are
    interpolated with the neighbouring background. The function identifies those affected edge
    voxels by comparing the difference of each voxel to its local neighbourhood. Background is
    assumed to be filled by zeroes and identified edge voxels are set to the background value.
    Inputs:
        *input_cmap: filename of 4d coordinate mapping.
        *path_output: path where output is saved.
        *basename_output: basename of output file.
        *edge_threshold: maximum difference to neighbouring voxel.
    
    created by Daniel Haenelt
    Date created: 26-10-2019
    Last modified: 27-10-2019
    """
    import os
    import numpy as np
    import nibabel as nb

    # load input
    cmap = nb.load(input_cmap)
    cmap_array = cmap.get_fdata()

    # get slab coordinates within volume (assumes background filled with zeroes)
    cmap_slab_coords = np.where(cmap_array != 0)

    # number of voxels
    n_vox = len(cmap_slab_coords[0])
    c_vox = 0
    
    # initialise for printing out loop status
    c_step = 0
    n_step = [10,20,30,40,50,60,70,80,90,100]

    # identify edges
    for n in range(n_vox):

        # current voxel
        i = cmap_slab_coords[0][n]
        j = cmap_slab_coords[1][n]
        k = cmap_slab_coords[2][n]
        t = cmap_slab_coords[3][n]
        
        # initialise voxel and its neighbours
        cmap_point = cmap_array[i,j,k,t]
        cmap_neighbour = []
            
        # consider neighour if not at edge of image volume
        if i > 0:
            cmap_neighbour.append(cmap_array[i-1,j,k,t])
        
        if j > 0:
            cmap_neighbour.append(cmap_array[i,j-1,k,t])
            
        if k > 0:
            cmap_neighbour.append(cmap_array[i,j,k-1,t])
            
        if i < np.shape(cmap_array)[0] - 1:
            cmap_neighbour.append(cmap_array[i+1,j,k,t])
    
        if j < np.shape(cmap_array)[1] - 1:
            cmap_neighbour.append(cmap_array[i,j+1,k,t])
    
        if k < np.shape(cmap_array)[2] - 1:
            cmap_neighbour.append(cmap_array[i,j,k+1,t])
    
        if i > 0 and j > 0:
            cmap_neighbour.append(cmap_array[i-1,j-1,k,t])
    
        if i > 0 and j < np.shape(cmap_array)[1] - 1:
            cmap_neighbour.append(cmap_array[i-1,j+1,k,t])
    
        if i > np.shape(cmap_array)[0] - 1 and j > 0:
            cmap_neighbour.append(cmap_array[i+1,j-1,k,t])
    
        if i > np.shape(cmap_array)[0] - 1 and j < np.shape(cmap_array)[1] - 1:
            cmap_neighbour.append(cmap_array[i+1,j+1,k,t])

        if i > 0 and k > 0:
            cmap_neighbour.append(cmap_array[i-1,j,k-1,t])
    
        if i > 0 and k < np.shape(cmap_array)[2] - 1:
            cmap_neighbour.append(cmap_array[i-1,j,k+1,t])
    
        if i > np.shape(cmap_array)[0] - 1 and k > 0:
            cmap_neighbour.append(cmap_array[i+1,j,k-1,t])
    
        if i > np.shape(cmap_array)[0] - 1 and j < np.shape(cmap_array)[2] - 1:
            cmap_neighbour.append(cmap_array[i+1,j,k+1,t])   

        if j > 0 and k > 0:
            cmap_neighbour.append(cmap_array[i,j-1,k-1,t])
    
        if j > 0 and k < np.shape(cmap_array)[2] - 1:
            cmap_neighbour.append(cmap_array[i,j-1,k+1,t])
    
        if j > np.shape(cmap_array)[1] - 1 and k > 0:
            cmap_neighbour.append(cmap_array[i,j+1,k-1,t])
    
        if j > np.shape(cmap_array)[1] - 1 and j < np.shape(cmap_array)[2] - 1:
            cmap_neighbour.append(cmap_array[i,j+1,k+1,t])   
             
        # compare voxel to its neighbours
        cmap_temp = np.abs(np.array(cmap_neighbour) - cmap_point)
        cmap_temp = cmap_temp[~np.isnan(cmap_temp)]
        cmap_temp = np.any(cmap_temp > edge_threshold)
            
        if cmap_temp:
            cmap_array[i,j,k,t] = np.nan
        
        # print status
        counter = np.floor(c_vox / n_vox * 100).astype(int)
        if counter == n_step[c_step]:
            print(basename_output+": "+str(counter)+" %")
            c_step += 1
        
        # counter loop cycles
        c_vox += 1

    # remove edges
    cmap_array[np.isnan(cmap_array)] = 0
    
    # get binary mask from single dimensions
    mask_array = cmap_array[:,:,:,0] * cmap_array[:,:,:,1] * cmap_array[:,:,:,2]
    mask_array[mask_array != 0] = 1

    # mask cmap
    cmap_array[:,:,:,0] = cmap_array[:,:,:,0] * mask_array
    cmap_array[:,:,:,1] = cmap_array[:,:,:,1] * mask_array
    cmap_array[:,:,:,2] = cmap_array[:,:,:,2] * mask_array

    # write output  
    output = nb.Nifti1Image(cmap_array, cmap.affine, cmap.header)
    if os.path.splitext(input_cmap)[1] == ".gz":
        nb.save(output,os.path.join(path_output,basename_output+"_edge.nii.gz"))
    else:
        nb.save(output,os.path.join(path_output,basename_output+"_edge.nii"))