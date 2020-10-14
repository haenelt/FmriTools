# -*- coding: utf-8 -*-

# python standard library inputs
import sys

# external inputs
import numpy as np
from numpy.fft import fft, ifft, fft2, ifft2

# local inputs
from fmri_tools.simulation.get_white import get_white_2d, get_white_1d
from fmri_tools.simulation.filter_bold import filter_bold_2d, filter_bold_1d
from fmri_tools.simulation.mask_pattern import mask_pattern_2d, mask_pattern_1d


def pattern_2d(Nx_sim=1024, Ny_sim=1024, FOVx=20, FOVy=20, Nx_mri=100, 
               Ny_mri=100, omega_x=0.5, omega_y=0.5, phi_x=0, phi_y=0, theta=0, 
               rect_shape=True, beta=0.05, fwhm_bold=1.02, fwhm_noise=0.001, 
               a_mask=1000, b_mask=1000, alpha_mask=0):
    """ Pattern 2D

    This function geometrical stripe pattern in 2D. The rest is similar to the 
    ODC pattern generation.    

    Parameters
    ----------
    Nx_sim : int, optional
        Array size of the simulated patch in x-direction. The default is 1024.
    Ny_sim : int, optional
        Array size of the simulated patch in y-direction. The default is 1024.
    FOVx : float, optional
        Field of view in x-direction (mm). The default is 20.
    FOVy : float, optional
        Field of view in y-direction (mm). The default is 20.
    Nx_mri : int, optional
        Array size of the MR image in x-direction. The default is 100.
    Ny_mri : int, optional
        Array size of the MR image in y-direction. The default is 100.
    omega_x : float, optional
        Frequency in x-direction in cycles/mm. The default is 0.5.
    omega_y : float, optional
        Frequency in y-direction in cycles/mm. The default is 0.5.
    phi_x : float, optional
        Phase in x-direction in deg. The default is 0.
    phi_y : float, optional
        Phase in y-direction in deg. The default is 0.
    theta : float, optional
        Rotation angle in deg. The default is 0.
    rect_shape : bool, optional
        Rectangular or sine pattern. The default is True.
    beta : float, optional
        Maximal BOLD response corresponding to neural response of 1. The default 
        is 0.05.
    fwhm_bold : float, optional
        BOLD point-spread width in mm. The default is 1.02.
    fwhm_noise : float, optional
        Measurement noise of BOLD response. The default is 0.001.
    a_mask : float, optional
        Major axis of elliptical mask. The default is 1000.
    b_mask : float, optional
        Minor axis of elliptical mask. The default is 1000.
    alpha_mask : float, optional
        Rotational angle of elliptical mask. The default is 0.

    Returns
    -------
    pattern_img : ndarray
        Neural map.
    y_img : ndarray
        BOLD response with measurement noise.
    ymri_img : ndarray
        Sampled MRI signal.
    F_bold_fft : ndarray
        BOLD filter in spatial frequency representation.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 07-01-2019     
    Last modified: 12-10-2020    

    """
    
    # generate geometrical pattern
    x = np.linspace(0,FOVx,Nx_sim)
    y = np.linspace(0,FOVy,Ny_sim)
    
    x_grid, y_grid = np.meshgrid(x,y)
    
    # compute rotation
    theta = np.radians(theta)
    c, s = np.cos(theta), np.sin(theta)
    x_rot = c*x_grid - s*y_grid
    y_rot = s*x_grid + c*y_grid
    
    # convert phase from deg to rad
    phi_x = np.radians(phi_x)
    phi_y = np.radians(phi_y)

    cos_x = np.cos(2*np.pi*omega_x*x_rot+phi_x)
    cos_y = np.cos(2*np.pi*omega_y*y_rot+phi_y)

    # merge pattern of x- and y-direction
    pattern_img = cos_x * cos_y

    # generate rect or sine pattern
    if rect_shape:
        pattern_img[pattern_img > 0] = 1
        pattern_img[pattern_img != 1] = -1

    # BOLD response
    F_bold_fft = filter_bold_2d(Nx_sim, Ny_sim, FOVx, FOVy, fwhm_bold, beta)
    y_img = np.real(ifft2( fft2(pattern_img) * F_bold_fft ))

    # add measurement noise
    noise_img = get_white_2d(Nx_sim, Ny_sim, 0, fwhm_noise/(2*np.sqrt(2*np.log(2))))
    y_img = noise_img + y_img

    # voxel sampling
    y_fft = fft2(y_img)

    kx_sample = np.round(Nx_mri/2).astype(int)
    ky_sample = np.round(Ny_mri/2).astype(int)

    ymri_fft = np.zeros((Nx_mri,Ny_mri), dtype=complex)
    ymri_fft[:kx_sample,:ky_sample] = y_fft[:kx_sample,:ky_sample] 
    ymri_fft[:kx_sample,-1:-ky_sample-1:-1] = y_fft[:kx_sample,-1:-ky_sample-1:-1] 
    ymri_fft[-1:-kx_sample-1:-1,:ky_sample] = y_fft[-1:-kx_sample-1:-1,:ky_sample] 
    ymri_fft[-1:-kx_sample-1:-1,-1:-ky_sample-1:-1] = y_fft[-1:-kx_sample-1:-1,-1:-ky_sample-1:-1] 

    ymri_img = np.real(ifft2(ymri_fft))
    
    # mask voxel sampling
    ymri_img = ymri_img * mask_pattern_2d(np.shape(ymri_img)[0], 
                                          np.shape(ymri_img)[1], 
                                          a_mask,
                                          b_mask,
                                          alpha_mask)
    
    return pattern_img, y_img, ymri_img, F_bold_fft


def pattern_1d(N_sim=1024, FOV=20, N_mri=100, omega=0.5, phi=0, rect_shape=True, 
               beta=0.05, fwhm_bold=1.02, fwhm_noise=0.001, a_mask=1000, 
               b_mask=1000):
    """ Pattern 1D
    
    This function geometrical stripe pattern in 1D. The rest is similar to the 
    stripe pattern generation in 2D.

    Parameters
    ----------
    N_sim : int, optional
        Array size of the simulated patch. The default is 1024.
    FOV : float, optional
        Field of view (mm). The default is 20.
    N_mri : int, optional
        Array size of the MR image. The default is 100.
    omega : float, optional
        Frequency in cycles/mm. The default is 0.5.
    phi : float, optional
        Phase in deg. The default is 0.
    rect_shape : bool, optional
        Rectangular or sine pattern. The default is True.
    beta : float, optional
        Maximal BOLD response corresponding to neural response of 1. The default 
        is 0.05.
    fwhm_bold : float, optional
        BOLD point-spread width in mm. The default is 1.02.
    fwhm_noise : float, optional
        Measurement noise of BOLD response. The default is 0.001.
    a_mask : float, optional
        Left side length of mask. The default is 1000.
    b_mask : float, optional
        Right side length of mask. The default is 1000.

    Returns
    -------
    pattern_img : ndarray
        Neural map.
    y_img : ndarray
        BOLD response with measurement noise.
    ymri_img : ndarray
        Sampled MRI signal.
    F_bold_fft : ndarray
        BOLD filter in spatial frequency representation.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 07-01-2019     
    Last modified: 12-10-2020
    
    """
    
    # add path of the executed script to the interpreter's search path
    sys.path.append(sys.argv[0])

    # generate geometrical pattern
    x = np.linspace(0,FOV,N_sim)
       
    # convert phase from deg to rad
    phi = np.radians(phi)

    pattern_img = np.cos(2*np.pi*omega*x+phi)

    # generate rect or sine pattern
    if rect_shape:
        pattern_img[pattern_img > 0] = 1
        pattern_img[pattern_img != 1] = -1

    # BOLD response
    F_bold_fft = filter_bold_1d(N_sim, FOV, fwhm_bold, beta)
    y_img = np.real(ifft( fft(pattern_img) * F_bold_fft ))

    # add measurement noise
    noise_img = get_white_1d(N_sim, 0, fwhm_noise/(2*np.sqrt(2*np.log(2))))
    y_img = noise_img + y_img

    # voxel sampling
    y_fft = fft(y_img)

    k_sample = np.round(N_mri/2).astype(int)

    ymri_fft = np.zeros(N_mri, dtype=complex)
    
    ymri_fft[:k_sample] = y_fft[:k_sample] 
    ymri_fft[-1:-k_sample-1:-1] = y_fft[-1:-k_sample-1:-1] 

    ymri_img = np.real(ifft(ymri_fft))
    
    # mask voxel sampling
    ymri_img = ymri_img * mask_pattern_1d(len(ymri_img), a_mask, b_mask)
    
    return pattern_img, y_img, ymri_img, F_bold_fft
