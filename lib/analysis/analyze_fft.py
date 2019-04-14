def analyze_fft(input, fovx, fovy, xv, yv, p_min=10, p_max=None, nsample=1000):
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
        *p_min: minimum prominence for peak detection.
        *p_max: maximum prominence for peak detection.
        *nsample: number of sampling point along projection line.
    Outputs:
        *P_max: maximum peak power relative to central frequency along projection axis.
        *k_max: corresponding spatial frequency.
        *k_line: spatial frequencies along projection axis.
        *fft_line: corresponding spectral power along projection axis.
        
    created by Daniel Haenelt
    Date created: 14-04-2019
    Last modified: 14-04-2019
    """
    import numpy as np
    from scipy.signal import find_peaks
    from lib.utils import get_fft

    # add one if nsample is an odd integer
    if np.mod(nsample,2) > 0:
        nsample += 1

    # get size of array
    x_size = np.shape(input)[0]
    y_size = np.shape(input)[1]
    
    # get coordinate of array
    kx = np.linspace(-1,1,x_size)*(x_size/(2*fovx))
    ky = np.linspace(-1,1,y_size)*(y_size/(2*fovy))
    kx_mesh, ky_mesh = np.meshgrid(kx,ky)

    # plot projection line from pca eigenvector. Because of nearest neighbor interpolation, the axis 
    # vector end point is <size>-1. The number of samples is hard coded.
    # get spatial frequency projection lines of major and minor axes from pca eigenvectors. The 
    # number of samples is hard coded.
    if np.abs(xv) < np.abs(yv):
        y_line = np.linspace(0,y_size-1,nsample)
        x_line = xv / yv * y_line + y_size/2
        x_line = x_line + x_size/2 - x_line[int(nsample/2)]
    else:
        x_line = np.linspace(0,x_size-1,nsample)
        y_line = yv / xv * x_line + x_size/2
        y_line = y_line + y_size/2 - y_line[int(nsample/2)]

    kxx_line = kx_mesh[np.round(y_line).astype(int),np.round(x_line).astype(int)]
    kyy_line = ky_mesh[np.round(y_line).astype(int),np.round(x_line).astype(int)]
   

    # get final k axes
    k_line = np.sqrt(kxx_line**2 + kyy_line**2)
    
    array_fft = get_fft(input)
    fft_line = array_fft[np.round(y_line).astype(int),np.round(x_line).astype(int)]
    
    # get one-sided spectrum
    arg_null = np.asscalar(np.argwhere(k_line==np.min(k_line))[-1])

    fft_line = fft_line[arg_null:]
    k_line = k_line[arg_null:]
       
    # normalize by central frequency
    fft_line = fft_line / fft_line[0] * 100
    
    # find peaks
    peak = find_peaks(fft_line, prominence=(p_min,p_max))[0]
    
    if len(peak) < 1:
        P_max = np.nan
        k_max = np.nan
    else:    
        P_max = np.asscalar(fft_line[peak[np.argwhere(fft_line[peak] == np.max(fft_line[peak]))]])
        k_max = np.asscalar(k_line[peak[np.argwhere(fft_line[peak] == np.max(fft_line[peak]))]])
    
    return k_max, P_max, k_line, fft_line