def analyze_fft(input, fovx, fovy, xv, yv, prominence_min=10, prominence_max=None):
    """
    This function computes the peak frequency and its corresponding power from a one-sided power 
    spectrum sampled on two projection lines (major and minor axes) with nearest neighbor 
    interpolation from a 2D array.
    Inputs:
        *input: input array (2D).
        *fovx: field of view in x-direction (mm).
        *fovy: field of view in y-firection (mm).
        *xv: x-coordinate of pca eigenvector.
        *yv: y-coordinate of pca eigenvector.
        *prominence_min: minimum prominence for peak detection.
        *prominence_max: maximum prominence for peak detection.
    Outputs:
        *P_max: maximum peak power relative to central frequency along projection axis.
        *k_max: corresponding spatial frequency.
        
    created by Daniel Haenelt
    Date created: 14.04.2019
    Last modified: 14.04.2019
    """
    import numpy as np
    from scipy.signal import find_peaks
    from lib.utils import get_fft

    # get size of array
    x_size = np.shape(input)[0]
    y_size = np.shape(input)[1]
    
    # get k-space representation of array
    kx = np.linspace(-x_size/2,x_size/2,x_size)/fovx
    ky = np.linspace(-y_size/2,y_size/2,y_size)/fovy
    kx_mesh, ky_mesh = np.meshgrid(kx,ky)    
    
    # plot projection line from pca eigenvector. Because of nearest neighbor interpolation, the axis 
    # vector end point is <size>-1. The number of samples is hard coded.
    # get spatial frequency projection lines of major and minor axes from pca eigenvectors. The 
    # number of samples is hard coded.
    if np.abs(xv) < np.abs(yv):
        y_line = np.linspace(0,y_size-1,1000)
        x_line = xv / yv * y_line + y_size/2
        x_line = x_line + x_size/2 - x_line[500]
        
        ky_line = np.linspace(-1,1,1000)*(y_size/fovy)
        kx_line = xv / yv * ky_line
    else:
        x_line = np.linspace(0,x_size-1,1000)
        y_line = yv / xv * x_line + x_size/2
        y_line = y_line + y_size/2 - y_line[500]

        kx_line = np.linspace(-1,1,1000)*(x_size/fovx)
        ky_line = yv / xv * kx_line
  
    # get final k axes
    k_line = np.sqrt(kx_line**2 + ky_line**2)
    
    array_fft = get_fft(input)
    fft_line = array_fft[np.round(y_line).astype(int),np.round(x_line).astype(int)]
    
    # get one-sided spectrum
    arg_null = np.asscalar(np.argwhere(k_line==np.min(k_line))[-1])

    fft_line = fft_line[arg_null:]
    k_line = k_line[arg_null:]
       
    # normalize by central frequency
    fft_line = fft_line / fft_line[0] * 100
    
    # find peaks
    peak = find_peaks(fft_line, prominence=(prominence_min,prominence_max))[0]
    
    if len(peak) < 1:
        P_max = np.nan
        k_max = np.nan
    else:    
        P_max = np.asscalar(fft_line[peak[np.argwhere(fft_line[peak] == np.max(fft_line[peak]))]])
        k_max = np.asscalar(k_line[peak[np.argwhere(fft_line[peak] == np.max(fft_line[peak]))]])
    
    return P_max, k_max