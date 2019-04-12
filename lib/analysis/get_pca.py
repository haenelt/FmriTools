def get_pca(input, n_iter=10, fft_threshold=0.25):
    """
    This function computes the principal axes from the thresholded fourier spectrum of an input 
    array. They are found by estimating the eigenvectors of the moments of inertia matrix I wich is
    defined by I = sum_i=0^N ( y_i^2 & -x_iy_i \\ -x_iy_i & x_i^2 ). x_i and y_i are the coordinates
    of the remaining points of the thresholded Fourier spectrum. This procedure mainly follows Borri 
    et al. (2016).
    Inputs:
        *input: input array.
        *n_iter: numer of iterations for fft normalization.
        *fft_threshold: thresholded for getting the coordinates of the central part.
    Outputs:
        *x_v1: x-coordinate of the major axis.
        *y_v1: y-coordinate of the major axis.
        *x_v2: x-coordinate of the minor axis.
        *y_v2: y-coordinate of the minor axis.
        
    created by Daniel Haenelt
    Date created: 12.04.2019
    Last modified: 12.04.2019
    """
    #import os
    import numpy as np
    from lib.utils import get_fft
    
    # compute normalized fourier spectrum of input array
    data_fft = get_fft(input, write_output = False, normalization = True, N = 10)
    
    # threshold spectrum
    data_fft[data_fft < fft_threshold] = 0
    data_fft[data_fft != 0] = 1

    # define mesh grid
    x = np.linspace(-np.shape(input)[0]/2,np.shape(input)[0]/2,np.shape(input)[0])
    y = np.linspace(-np.shape(input)[1]/2,np.shape(input)[1]/2,np.shape(input)[1])
    x_mesh, y_mesh = np.meshgrid(x,y)

    # get coordinates of thresholded spectrum
    x_mesh = x_mesh * data_fft
    y_mesh = y_mesh * data_fft

    x_mesh = x_mesh[x_mesh != 0]
    y_mesh = y_mesh[y_mesh != 0]

    # compute moments of inertia
    a = np.sum(x_mesh**2)
    b = np.sum(y_mesh**2)
    c = np.sum(x_mesh*y_mesh)

    I = np.array([[b, -c],[-c, a]])

    # compute eigenvalues and eigenvectors of matrix I
    evals, evecs = np.linalg.eig(I)

    # sort eigenvalues in decreasing order
    sort_indices = np.argsort(evals)[::-1]
    x_v1, y_v1 = evecs[:, sort_indices[0]] # Eigenvector with largest eigenvalue
    x_v2, y_v2 = evecs[:, sort_indices[1]]
    
    return x_v1, y_v1, x_v2, y_v2
