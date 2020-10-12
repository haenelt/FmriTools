# -*- coding: utf-8 -*-

# python standard library inputs
import sys

# external inputs
import numpy as np
import nibabel as nb
from numpy.fft import fft, ifft, fft2, ifft2

# local inputs
from fmri_tools.simulation.get_white import get_white_2d, get_white_1d
from fmri_tools.simulation.filter_odc import filter_odc_2d, filter_odc_1d
from fmri_tools.simulation.filter_bold import filter_bold_2d, filter_bold_1d
from fmri_tools.simulation.mask_pattern import mask_pattern_2d, mask_pattern_1d
from fmri_tools.simulation.filter_sigmoid import filter_sigmoid


def odc_2d(Nx_sim=1024, Ny_sim=1024, FOVx=20, FOVy=20, Nx_mri=100, Ny_mri=100, 
           rho=0.5, delta=0.3, epsilon=0.4, theta=0, alpha=4, beta=0.05, 
           fwhm_bold=1.02, fwhm_noise=0.001, a_mask=1000, b_mask=1000, 
           alpha_mask=0, path_white=False):
    """
    This function generates realistic ocular dominance patterns according to the 
    model proposed by Rojer and Schwartz (1990) in 2D. Implementation follows 
    (Chaimow et al., 2011) and most of the default values for the columnar 
    pattern are taken from this publication. First, a white noise pattern is 
    defined and an anisotropic band-pass filter and a non-linear sigmoidal 
    filter are applied. This neural map is converted to a BOLD response map by 
    applying a Gaussian filter. White noise is added to simulate uncorrelated 
    measurement noise to the signal. The MRI sampling procedure is considered by 
    only inner k-space lines of the spatial frequency representation of the BOLD 
    reseponse map. 
    Inputs:
        *Nx_sim: array size of the simulated patch in x-direction (use only even integers).
        *Ny_sim: array size of the simulated patch in y-direction (use only even integers).
        *FOVx: field of view in x-direction (mm).
        *FOVy: field of view in y-direction (mm).
        *Nx_mri: array size of the MR image in x-direction.
        *Ny_mri: array size of the MR image in y-direction.
        *rho: main spatial frequency determining columnar width in cycles/mm.
        *delta: variations orthogonal to ODC bands in cycles/mm (irregularity).
        *epsilon: variations parallel to ODC bands in cycles/mm (branchiness).
        *theta: orientation of the columnar pattern in deg.
        *alpha: sharpness parameter of the sigmoidal filter.
        *beta: maximal BOLD response corresponding to neural response of 1.
        *fwhm_bold: BOLD point-spread width in mm.
        *fwhm_noise: measurement noise of BOLD response.
        *a_mask: major axis of elliptical mask.
        *b_mask: minor axis of elliptical mask.
        *alpha_mask: rotational angle of elliptical mask.
        *path_white: path to existing white noise image.
    Outputs:
        *white_img: initial white noise pattern.
        *odc_img: neural map.
        *y_img: BOLD response with measurement noise.
        *ymri_img = sampled MRI signal.
        *F_odc_fft: anisotropic band-pass filter in spatial frequency representation.
        *F_bold_fft: BOLD filter in spatial frequency representation.
        
    created by Daniel Haenelt
    Date created: 07-01-2019     
    Last modified: 12-10-2020
    """    

    # get white noise
    if path_white:
        input = nb.load(path_white)
        white_img = input.get_fdata()
    else:     
        white_img = get_white_2d(Nx_sim, Ny_sim, 0, 1)

    # get band-pass filter
    F_odc_fft = filter_odc_2d(Nx_sim, Ny_sim, FOVx, FOVy, rho, delta, epsilon, theta)

    # generate ODC pattern (neural map)
    odc_img = np.real(ifft2( fft2(white_img) * F_odc_fft ))
    odc_img = filter_sigmoid(odc_img, alpha)

    # BOLD response
    F_bold_fft = filter_bold_2d(Nx_sim, Ny_sim, FOVx, FOVy, fwhm_bold, beta)
    y_img = np.real(ifft2( fft2(odc_img) * F_bold_fft ))

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
    
    return white_img, odc_img, y_img, ymri_img, F_odc_fft, F_bold_fft


def odc_1d(N_sim=1024, FOV=20, N_mri=100, rho=0.5, delta=0.3, alpha=4, 
           beta=0.05, fwhm_bold=1.02, fwhm_noise=0.001, a_mask=1000, 
           b_mask=1000, path_white=False):
    """
    This function generates realistic ocular dominance patterns according to the 
    model proposed by Rojer and Schwartz (1990) in 1D. The rest follows similar 
    to the ODC generation in 2D.
    Inputs:
        *N_sim: array size of the simulated patch in x-direction (use only even integers).
        *FOV: field of view in x-direction (mm).
        *N_mri: array size of the MR image in x-direction.
        *rho: main spatial frequency determining columnar width in cycles/mm.
        *delta: variations orthogonal to ODC bands in cycles/mm (irregularity).
        *alpha: sharpness parameter of the sigmoidal filter.
        *beta: maximal BOLD response corresponding to neural response of 1.
        *fwhm_bold: BOLD point-spread width in mm.
        *fwhm_noise: measurement noise of BOLD response.
        *a_mask: left side length of mask.
        *b_mask: right side length of mask.
        *path_white: path to existing white noise image.
    Outputs:
        *white_img: initial white noise pattern.
        *odc_img: neural map.
        *y_img: BOLD response with measurement noise.
        *ymri_img = sampled MRI signal.
        *F_odc_fft: anisotropic band-pass filter in spatial frequency representation.
        *F_bold_fft: BOLD filter in spatial frequency representation.
        
    created by Daniel Haenelt
    Date created: 08-01-2019     
    Last modified: 12-10-2020
    """
    
    # add path of the executed script to the interpreter's search path
    sys.path.append(sys.argv[0])
    
    # get white noise
    if path_white:
        input = nb.load(path_white)
        white_img = input.get_fdata()
    else:     
        white_img = get_white_1d(N_sim, 0, 1)

    # get band-pass filter
    F_odc_fft = filter_odc_1d(N_sim, FOV, rho, delta)

    # generate ODC pattern (neural map)
    odc_img = np.real(ifft( fft(white_img) * F_odc_fft ))
    odc_img = filter_sigmoid(odc_img, alpha)

    # BOLD response
    F_bold_fft = filter_bold_1d(N_sim, FOV, fwhm_bold, beta)
    y_img = np.real(ifft( fft(odc_img) * F_bold_fft ))

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
    
    return white_img, odc_img, y_img, ymri_img, F_odc_fft, F_bold_fft
