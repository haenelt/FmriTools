def get_bandpass_filter(nx, ny, fovx, fovy, kcut_low=0, kcut_high=1, apply_fftshift=False, k_fwhm=0, 
                        theta_fwhm=0, theta1=0, theta2=180):
    """
    This function generates a bandpass filter for an input image defined by lower and upper cutoff
    frequencies. Additionally, Gaussian attenuation of border can be appled and the filter can be 
    only defined to a specific k-space region defined by lower and upper angle thresholds. 
    Inputs:
        *nx: matrix size of input array in x-direction.
        *ny: matrix size of input array in x-direction.
        *fovx: field of view of input array in x-direction in mm.
        *fovy: field of view of input array in y-direction in mm.
        *kcut_low: lower spatial cutoff frequency in cycles/mm.
        *kcut_high: upper spatial cutoff frequency in cycles/mm.
        *apply_fftshift: fftshift to stay in the convention of spatial frequencies in numpy arrays.
        *k_fwhm: full-width at half maximum of gaussian filter in frequency direction.
        *theta_fwhm: full-width at half maximum of gaussian filter in angle direction.
        *theta1: lower cutoff angle in deg [0,180].
        *theta2: higher cutoff angle in deg [0,180].
    Outputs:
        *B: spatial frequency filter.
        
    created by Daniel Haenelt
    Date created: 31.07.2019         
    Last modified: 31.07.2019
    """
    import numpy as np
    from numpy.fft import fftshift

    # parameters of gaussian
    beta = 1
    k_sigma = k_fwhm / ( 2*np.sqrt(2*np.log(2)) )
    theta_sigma = theta_fwhm / ( 2*np.sqrt(2*np.log(2)) )

    # get maximum k-space coordinate in x- and y-direction
    kx_max = nx/(2*fovx)
    ky_max = ny/(2*fovy)
    
    # get k-space axes
    kx = np.linspace(-kx_max,kx_max,nx)
    ky = np.linspace(-ky_max,ky_max,ny)
    
    # define two-dimensional k-space grid
    kx_grid, ky_grid = np.meshgrid(kx,ky)
    
    # convert to polar coordinates
    K_r = np.sqrt(kx_grid**2 + ky_grid**2)
    
    K_pol = np.arctan2(ky_grid,kx_grid)
    K_pol[K_pol < 0] = K_pol[K_pol < 0] + np.pi
    K_pol = K_pol / np.pi * 180
    
    # frequency filter
    B_r = K_r.copy()
    B_r[np.where(np.logical_and(B_r>=kcut_low, B_r<=kcut_high))] = np.nan
    B_r[~np.isnan(B_r)] = 0
    B_r[B_r != 0] = 1
    
    # angle filter
    B_pol = K_pol.copy()
    if theta2 > theta1:
        B_pol[np.where(np.logical_and(B_pol>=theta1, B_pol<=theta2))] = np.nan
        B_pol[~np.isnan(B_pol)] = 0
        B_pol[B_pol != 0] = 1
    else:
        B_pol[np.where(np.logical_and(B_pol<=theta1, B_pol>=theta2))] = np.nan
        B_pol[np.isnan(B_pol)] = 0
        B_pol[B_pol != 0] = 1
        B_pol[K_pol == 0] = 1 # important to also fill the phase gap
    
    # filter edges
    if k_fwhm != 0:        
        B1 = beta / (np.sqrt(2*np.pi)*k_sigma) * np.exp( -(K_r-kcut_low)**2 / ( 2 * k_sigma**2 ) )
        B2 = beta / (np.sqrt(2*np.pi)*k_sigma) * np.exp( -(K_r-kcut_high)**2 / ( 2 * k_sigma**2 ) )
        
        B_temp = B1 + B2
        B_temp2 = B_temp.copy()
        B_temp2[B_r == 1] = 0
        B_max = np.max(B_temp2) 
        B_temp[B_r == 1] = B_max
        B_r = B_temp.copy()
        
        # normalize filter
        B_r = ( B_r - np.min(B_r) ) / ( np.max(B_r) - np.min(B_r) )   
    
    if theta_fwhm != 0 and theta2-theta1 < 180:
        angle1 = np.mod(K_pol-theta1,180)
        angle2 = np.mod(theta1-K_pol,180)
        angle3 = np.mod(K_pol-theta2,180)
        angle4 = np.mod(theta2-K_pol,180)        
        
        B_pol1a = beta / (np.sqrt(2*np.pi)*theta_sigma) * np.exp( -angle1**2 / ( 2 * theta_sigma**2 ) )
        B_pol1b = beta / (np.sqrt(2*np.pi)*theta_sigma) * np.exp( -angle2**2 / ( 2 * theta_sigma**2 ) )
        B_pol2a = beta / (np.sqrt(2*np.pi)*theta_sigma) * np.exp( -angle3**2 / ( 2 * theta_sigma**2 ) )
        B_pol2b = beta / (np.sqrt(2*np.pi)*theta_sigma) * np.exp( -angle4**2 / ( 2 * theta_sigma**2 ) )

        B_temp = B_pol1a + B_pol1b + B_pol2a + B_pol2b
        B_temp2 = B_temp.copy()
        B_temp2[B_pol == 1] = 0
        B_max = np.max(B_temp2)
        B_temp[B_pol == 1] = B_max
        B_pol = B_temp.copy()
        
        # normalize filter
        B_pol = ( B_pol - np.min(B_pol) ) / ( np.max(B_pol) - np.min(B_pol) )

    # compose frequency and angle filter
    B = B_r * B_pol

    # shift zero k-space line to the left border to match the numpy fft convention
    if apply_fftshift == True:
        B = fftshift(B)
    
    return B