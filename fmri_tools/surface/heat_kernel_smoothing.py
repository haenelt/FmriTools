# -*- coding: utf-8 -*-

# external inputs
import numpy as np
from gbb.neighbor.nn_2d import nn_2d


def heat_kernel_smoothing(vtx, fac, data, adjm, sigma, n_smooth):
    """ Heat kernel smoothing
    
    This function performs heat kernel smoothing [1,2,3] on a triangle mesh. The 
    code is mainly adapted from the matlab code by Chung et al. [4]. The kernel 
    bandwidth corresponds to diffusion time in the heat equation [3]. The FWHM 
    follows 4*sqrt(log 2*n_smooth*sigma) with the natural log.
    
    If you use this code, please reference one of the following papers. The 
    details on the mathematical basis of of the algorithm can be found in these 
    papers.

    Parameters
    ----------
    vtx : ndarray
        Vertex points of surface mesh.
    fac : ndarray
        Faces of surface mesh.
    data : ndarray
        Array of vertex-wise sampled data points.
    adjm : ndarray
        Adjacency matrix.
    sigma : float
        Kernel bandwidth.
    n_smooth : int
        Number of iterations.

    Returns
    -------
    res : ndarray
        Array of vertex-wise smoothed data points.

    References
    -------
    .. [1] Chung, MK, et al. Cortical thickness analysis in autism via heat 
    kernel smoothing. Neuroimage 25(1), 1256--1265 (2005).
    .. [2] Chung, MK, et al. Unified statistical approach to cortical thickness 
    analysis. Inf Process Med Imaging 19, 627--638 (2005).
    .. [3] Chung, MK, et al. Encoding cortical surface by spherical harmonics. 
    Statistica Sinica 18, 1269--1291 (2008).
    .. [4] http://pages.stat.wisc.edu/~mchung/softwares/hk/hk_smooth.m

    Notes
    -------
    created by Daniel Haenelt
    Date created: 04-03-2020
    Last modified: 12-10-2020
    
    """
    
    # number of vertices
    n_vertex = len(vtx)

    # heat kernel shape
    K = lambda x, sigma : np.exp(-x/(4*sigma))/np.sum(np.exp(-x/(4*sigma)))
    
    # get max degree (number of first order neighbors)
    max_degree = 0
    for i in range(n_vertex):
        nn = nn_2d(i, adjm, 0)
        degree = len(nn)
        if degree > max_degree:
            max_degree = degree
    
    # heat kernel weight computation
    neighbor = np.zeros((n_vertex, max_degree+1)).astype(int) # including the current vertex
    weight = np.zeros((n_vertex, max_degree+1)) # including the current vertex
    for i in range(n_vertex):
        
        # get vertex neighbors
        nn = nn_2d(i, adjm, 0)
        degree = len(nn)

        # get distance to vertex neighbors
        distance = 0
        for j in range(degree):
            distance = np.append(distance, np.sum(( vtx[nn[j]] - vtx[i,:] ) ** 2))
            
        # get heat kernel weighting for each neighbor
        weight[i,:1+degree] = K(distance, sigma)
        
        # get corresponding neighbor (add 1 because of dummy row)
        neighbor[i,:1+degree] = np.append([i],nn) + 1        

    # add dummy row
    data = np.append(1, data)
    
    # iterative kernel smoothing
    for i in range(n_smooth):
    
        # add weights
        res = np.zeros_like(data)
        for j in range(max_degree):
            res[1:] += data[neighbor[:,j]] * weight[:,j]
        
        # initialize new data array
        data = res.copy()
    
    # remove dummy row    
    res = res[1:]
    
    return res
