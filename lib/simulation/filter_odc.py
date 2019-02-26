def filter_odc_2d(nx, ny, fovx, fovy, rho, delta, epsilon, theta):
    """
    This function defines the 2D ocular dominance column filter (anisotropic band-pass filter) in 
    spatial frequency space which is taken from Chaimow et al. (2016a). The filter is normalised in
    order to have the same variance as the white noise source image.
    Input:
        *nx: array size in x-direction.
        *ny: array size in y-direction.
        *fovx: field of view in x-direction (mm).
        *fovy: field of view in y-direction (mm).
        *rho: main spatial frequency determining columnar width in cycles/mm.
        *delta: variations orthogonal to ODC bands (irregularity).
        *epsilon: variations parallel to ODC bands (branchiness).
        *theta: rotation angle in deg.
    Outputs:
        * F: band-pass filter array in spatial frequency space.
        
    created by Daniel Haenelt
    Date created: 07-01-2019
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
    
    # convert theta from deg to rad
    theta = np.radians(theta)
    
    # convert to polar coordinates
    k_r = np.sqrt(kx_grid**2 + ky_grid**2)
    k_pol = np.arctan2(ky_grid,kx_grid)

    # ODC filter
    F_rad1 = np.exp( -(k_r-rho)**2 / (2*delta**2) )
    F_rad2 = np.exp( -(k_r+rho)**2 / (2*delta**2) )
    
    F_ang1 = np.exp( np.cos(k_pol-theta) / epsilon**2 )
    F_ang2 = np.exp( -np.cos(k_pol-theta) / epsilon**2 )

    F = (F_rad1 + F_rad2) * (F_ang1 + F_ang2)
    
    # normalise filter    
    c = np.sqrt( np.sum(F**2) / np.size(F) )
    
    F = F / c  
    
    # shift zero k-space line to the left border to match the numpy fft convention
    F = fftshift(F)
    
    return F

def filter_odc_1d(n, fov, rho, delta):
    """
    This function defines the 1D ocular dominance column filter (anisotropic band-pass filter) in 
    spatial frequency space which is taken from Chaimow et al. (2016a). The filter is normalised in
    order to have the same variance as the white noise source image.
    Input:
        *n: array size.
        *fov: field of view (mm).
        *rho: main spatial frequency determining columnar width in cycles/mm.
        *delta: variations orthogonal to ODC bands (irregularity).
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
    
    # ODC filter
    F_rad1 = np.exp( -(k-rho)**2 / (2*delta**2) )
    F_rad2 = np.exp( -(k+rho)**2 / (2*delta**2) )

    F = F_rad1 + F_rad2
    
    # normalise filter    
    c = np.sqrt( np.sum(F**2) / np.size(F) )
    
    F = F / c  
    
    # shift zero k-space line to the left border to match the numpy fft convention
    F = fftshift(F)
    
    return F