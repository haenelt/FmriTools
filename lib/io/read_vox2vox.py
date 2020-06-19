def read_vox2vox(input_lta):
    """
    This function reads a freesurfer lta file and extracts the vox2vox transformation matrix as 
    numpy array.
    Inputs:
        *input_lta: freesurfer lta file.
    Outputs:
        *M: forwards affine transformation matrix.
        *Minv: inverse of M.
            
    created by Daniel Haenelt
    Date created: 19-06-2020
    Last modified: 19-06-2020
    """
    import numpy as np

    with open(input_lta, "r") as f:
        x = f.readlines()

        # it is assumed that the vox2vox transformation matrix is found at specific lines in the 
        # lta file
        transformation = x[8:12]

    # convert matrix to numpy array
    M = np.zeros((4,4))
    for i in range(len(transformation)):
        x = transformation[i].split()
    
        M[i,0] = float(x[0])
        M[i,1] = float(x[1])
        M[i,2] = float(x[2])
        M[i,3] = float(x[3])
        
        # get inverted transformation matrix (pseudoinverse)
        Minv = np.linalg.pinv(M)
    
    return M, Minv