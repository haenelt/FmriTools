def pattern_2d(Nx_sim=1024, Ny_sim=1024, FOVx=20, FOVy=20, Nx_mri=100, Ny_mri=100, omega_x=0.5,
               omega_y=0.5, phi_x=0, phi_y=0, theta=0, rect_shape=True, beta=0.05, fwhm_bold=1.02, 
               fwhm_noise=0.001, a_mask=1000, b_mask=1000, alpha_mask=0):
    """
    This function geometrical stripe pattern in 2D. The rest is similar to the ODC pattern 
    generation.
    Inputs:
        *Nx_sim: array size of the simulated patch in x-direction (use only even integers).
        *Ny_sim: array size of the simulated patch in y-direction (use only even integers).
        *FOVx: field of view in x-direction (mm).
        *FOVy: field of view in y-direction (mm).
        *Nx_mri: array size of the MR image in x-direction.
        *Ny_mri: array size of the MR image in y-direction.
        *omega_x: frequency in x-direction in cycles/mm.
        *omega_y: frequency in y-direction inc cycles/mm.
        *phi_x: phase in x-direction in deg.
        *phi_y: phase in y-direction in deg.
        *theta: rotation angle in deg.
        *rect_shape (boolean): rectangular or sine pattern.
        *beta: maximal BOLD response corresponding to neural response of 1.
        *fwhm_bold: BOLD point-spread width in mm.
        *fwhm_noise: measurement noise of BOLD response.
        *a_mask: major axis of elliptical mask.
        *b_mask: minor axis of elliptical mask.
        *alpha_mask: rotational angle of elliptical mask.
    Outputs:
        *pattern_img: neural map.
        *y_img: BOLD response with measurement noise.
        *ymri_img: sampled MRI signal.
        *F_bold_fft: BOLD filter in spatial frequency representation.
        
    created by Daniel Haenelt
    Date created: 07-01-2019     
    Last modified: 07-01-2019
    """    
    import numpy as np
    from numpy.fft import fft2, ifft2
    from lib.simulation.get_white import get_white_2d
    from lib.simulation.filter_bold import filter_bold_2d
    from lib.simulation.mask_pattern import mask_pattern_2d

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

def pattern_1d(N_sim=1024, FOV=20, N_mri=100, omega=0.5, phi=0, rect_shape=True, beta=0.05, 
               fwhm_bold=1.02, fwhm_noise=0.001, a_mask=1000, b_mask=1000):
    """
    This function geometrical stripe pattern in 1D. The rest is similar to the stripe pattern 
    generation in 2D.
    Inputs:
        *N_sim: array size of the simulated patch (use only even integers).
        *FOV: field of view (mm).
        *N_mri: array size of the MR image.
        *omega: frequency in cycles/mm.
        *phi: phase in deg.
        *theta: rotation angle in deg.
        *rect_shape (boolean): rectangular or sine pattern.
        *beta: maximal BOLD response corresponding to neural response of 1.
        *fwhm_bold: BOLD point-spread width in mm.
        *fwhm_noise: measurement noise of BOLD response.
        *a_mask: left side length of mask.
        *b_mask: right side length of mask.
    Outputs:
        *pattern_img: neural map.
        *y_img: BOLD response with measurement noise.
        *ymri_img: sampled MRI signal.
        *F_bold_fft: BOLD filter in spatial frequency representation.
        
    created by Daniel Haenelt
    Date created: 07-01-2019     
    Last modified: 07-01-2019
    """   
    # add path of the executed script to the interpreter's search path
    import sys
    sys.path.append(sys.argv[0])
    
    import numpy as np
    from numpy.fft import fft, ifft
    from lib.simulation.get_white import get_white_1d
    from lib.simulation.filter_bold import filter_bold_1d
    from lib.simulation.mask_pattern import mask_pattern_1d

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
