def filter_bold_2d(nx, ny, fovx, fovy, fwhm, beta):
    """
    This defines the BOLD blurring in 2D. It is defined as modulation transfer function in spatial 
    frequency space. The filter is the Fourier transform of a two-dimensional Gaussian kernel in 
    image space.
    Inputs:
        *nx: array size in x-direction.
        *ny: array size in y-direction.
        *fovx: field of view in x-direction (mm).
        *fovy: field of view in y-direction (mm).
        *fwhm: full width at half maximum of the Gaussian kernel.
        *beta: maximum BOLD response for neural response of 1.
    Outputs:
        * F: band-pass filter array in spatial frequency space.
        
    created by Daniel Haenelt
    Date created: 04-01-2019          
    Last modified: 07-01-2019
    """
    import numpy as np
    from numpy.fft import fftshift 

    # get maximum k-space coordinate in x- and y-direction
    kx_max = nx/(2*fovx)
    ky_max = ny/(2*fovy)
    
    # get k-space axes
    kx = np.linspace(-kx_max,kx_max,nx)
    ky = np.linspace(-ky_max,ky_max,ny) 
    
    # define two-dimensional k-space grid
    kx_grid, ky_grid = np.meshgrid(kx,ky)
    
    # convert to polar coordinates
    k_r = np.sqrt(kx_grid**2 + ky_grid**2)
       
    # BOLD modulation transfer function
    sigma = fwhm / ( 2*np.sqrt(2*np.log(2)) )
    
    F = beta * np.exp( -2 * np.pi**2 * sigma**2 * k_r**2 )

    # shift zero k-space line to the left border to match the numpy fft convention
    F = fftshift(F)
    
    return F

def filter_bold_1d(n, fov, fwhm, beta):
    """
    This defines the BOLD blurring in 1D. It is defined as modulation transfer function in spatial 
    frequency space. The filter is the Fourier transform of a one-dimensional Gaussian kernel in 
    image space.
    Inputs:
        *n: array size.
        *fov: field of view (mm).
        *fwhm: full width at half maximum of the Gaussian kernel.
        *beta: maximum BOLD response for neural response of 1.
    Outputs:
        * F: band-pass filter array in spatial frequency space.
        
    created by Daniel Haenelt
    Date created: 08-01-2019          
    Last modified: 08-01-2019
    """
    import numpy as np
    from numpy.fft import fftshift 

    # get maximum k-space coordinate in x- and y-direction
    k_max = n/(2*fov)
    
    # get k-space axes
    k = np.linspace(-k_max,k_max,n) 
       
    # BOLD modulation transfer function
    sigma = fwhm / ( 2*np.sqrt(2*np.log(2)) )
    
    F = beta * np.exp( -2 * np.pi**2 * sigma**2 * k**2 )

    # shift zero k-space line to the left border to match the numpy fft convention
    F = fftshift(F)
    
    return F