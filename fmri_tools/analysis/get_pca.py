# -*- coding: utf-8 -*-

# external inputs
import numpy as np

# local inputs
from ..utils import get_fft


def get_pca(arr, fft_threshold=0.25):
    """Get  PCA.
    
    This function computes the principal axes from the thresholded fourier 
    spectrum of an input array. They are found by estimating the eigenvectors of 
    the moments of inertia matrix I wich is defined by I = sum_i=0^N ( y_i^2 & 
    -x_iy_i \\ -x_iy_i & x_i^2 ). x_i and y_i are the coordinates of the 
    remaining points of the thresholded Fourier spectrum. This procedure mainly 
    follows [1].

    Parameters
    ----------
    arr : ndarray
        Input array.
    fft_threshold : float, optional
        Thresholded for getting the coordinates of the central part. The 
        default is 0.25.

    Returns
    -------
    x_v1 : float
        x-coordinate of the major axis.
    y_v1 : float
        y-coordinate of the major axis.
    x_v2 : float
        x-coordinate of the minor axis.
    y_v2 : float
        y-coordinate of the minor axis.
    
    References
    -------
    .. [1] Borri, Marco, et al. A novel approach to evaluate spatial resolution 
    of MRI clinical images for optimization and standardization of breast 
    screening protocols, Med Phys 43(12), 6354--6363 (2016).
    
    """

    # compute normalized fourier spectrum of input array
    data_fft = get_fft(arr, write_output=False, normalization=True, n=10)

    # threshold spectrum
    data_fft[data_fft < fft_threshold] = 0
    data_fft[data_fft != 0] = 1

    # define mesh grid
    x = np.linspace(-np.shape(arr)[0] / 2, np.shape(arr)[0] / 2,
                    np.shape(arr)[0])
    y = np.linspace(-np.shape(arr)[1] / 2, np.shape(arr)[1] / 2,
                    np.shape(arr)[1])
    y_mesh, x_mesh = np.meshgrid(y, x)

    # get coordinates of thresholded spectrum    
    x_mesh = x_mesh[data_fft != 0]
    y_mesh = y_mesh[data_fft != 0]

    # compute moments of inertia
    a = np.sum(x_mesh ** 2)
    b = np.sum(y_mesh ** 2)
    c = np.sum(x_mesh * y_mesh)

    m = np.array([[b, -c], [-c, a]])

    # compute eigenvalues and eigenvectors of matrix I
    evals, evecs = np.linalg.eig(m)

    # sort eigenvalues in decreasing order
    sort_indices = np.argsort(evals)[::-1]
    x_v1, y_v1 = evecs[:, sort_indices[0]]  # Eigenvector with largest eigenvalue
    x_v2, y_v2 = evecs[:, sort_indices[1]]

    return x_v1, y_v1, x_v2, y_v2
