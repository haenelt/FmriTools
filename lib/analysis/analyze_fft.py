def analyze_fft(input, fovx, fovy, x_v1, y_v1, x_v2, y_v2):
    """
    This function computes the peak frequency from a one-sided power spectrum taken from a 
    projection line of a 2D FFT array.
    Inputs:
        *input: input array (2D).
    Outputs:
        * bla: bla
        
    created by Daniel Haenelt
    Date created:           
    Last modified: 
    """
    import numpy as np
    from scipy.signal import find_peaks
    from lib.utils import get_fft

    prominence_min = 10
    prominence_max = None


    # get size of array
    x_size = np.shape(input)[0]
    y_size = np.shape(input)[1]
    
    # get k-space representation of array
    kx = np.linspace(-x_size/2,x_size/2,x_size)/fovx
    ky = np.linspace(-y_size/2,y_size/2,y_size)/fovy
    kx_mesh, ky_mesh = np.meshgrid(kx,ky)    
    
    # plot projection lines of major and minor axes from pca eigenvectors
    x_major = np.linspace(0,x_size-1,1000)
    y_major = y_v1 / x_v1 * x_major + x_size/2
    y_minor = np.linspace(0,y_size-1,1000)
    x_minor = x_v2 / y_v2 * y_minor + y_size/2
    
    # get minor and major axes from pca eigenvectors
    kx_major = kx_mesh[np.round(x_major).astype(int),np.round(y_minor).astype(int)]
    ky_major = ky_mesh[np.round(x_major).astype(int),np.round(y_major).astype(int)]
    kx_minor = kx_mesh[np.round(x_minor).astype(int),np.round(y_minor).astype(int)]
    ky_minor = ky_mesh[np.round(x_minor).astype(int),np.round(y_minor).astype(int)]
    
    # get final k axes
    k_major = np.sqrt(kx_major**2 + ky_major**2)
    k_minor = np.sqrt(kx_minor**2 + ky_minor**2)
    
    array_fft = get_fft(input)
    fft_major = array_fft[np.round(y_major).astype(int),np.round(x_major).astype(int)]
    fft_minor = array_fft[np.round(y_minor).astype(int),np.round(x_minor).astype(int)]
    
    # get one-sided spectrum
    major_arg_null = np.asscalar(np.argwhere(k_major==np.min(k_major))[-1])
    minor_arg_null = np.asscalar(np.argwhere(k_minor==np.min(k_minor))[-1])

    fft_major = fft_major[major_arg_null:]
    fft_minor = fft_minor[minor_arg_null:]
    k_major = k_major[major_arg_null:]
    k_minor = k_minor[minor_arg_null:]
    
    # normalize by central frequency
    fft_major = fft_major / fft_major[0] * 100
    fft_minor = fft_minor / fft_minor[0] * 100
    
    # find peaks
    major_peak = find_peaks(fft_major, prominence=(prominence_min,prominence_max))[0]
    minor_peak = find_peaks(fft_minor, prominence=(prominence_min,prominence_max))[0]
    
    if len(major_peak) < 1:
        P_max_major = np.nan
        k_max_major = np.nan
    else:    
        P_max_major = np.asscalar(fft_major[major_peak[np.argwhere(fft_major[major_peak] == np.max(fft_major[major_peak]))]])
        k_max_major = np.asscalar(k_major[major_peak[np.argwhere(fft_major[major_peak] == np.max(fft_major[major_peak]))]])
    
    if len(minor_peak) < 1:
        P_max_minor = np.nan
        k_max_minor = np.nan
    else:    
        P_max_minor = np.asscalar(fft_minor[minor_peak[np.argwhere(fft_minor[minor_peak] == np.max(fft_minor[minor_peak]))]])
        k_max_minor = np.asscalar(k_minor[minor_peak[np.argwhere(fft_minor[minor_peak] == np.max(fft_minor[minor_peak]))]])
    
    return P_max_major, k_max_major, P_max_minor, k_max_minor