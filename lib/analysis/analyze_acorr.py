def analyze_acorr(input, fovx, fovy, xv, yv, p_min=0.01, p_max=None, nsample=1000):
    """
    This function blas computes the normalized autocorrelation (NAC) from a 2D input array and
    estimates the width of the central peak and the distance to its first neighbor peak along two
    projection lines sampled with nearest neighbor interpolation.
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
        *fwhm_central: fwhm of nac central peak in mm along projection axis.
        *P_neighbor: power of first neighbor peak along projection axis.
        *d_neighbor: distance to first neighbor peak in mm along projection axis.
        *d: lag in mm along projection axis.
        *acorr_line: nac long projection axis.
        
    created by Daniel Haenelt
    Date created: 14-04-2019
    Last modified: 14-04-2019
    """
    import numpy as np
    from scipy.signal import find_peaks
    from lib.utils import get_acorr

    # add one if nsample is an odd integer
    if np.mod(nsample,2) > 0:
        nsample += 1

    # get size of array
    x_size = np.shape(input)[0]
    y_size = np.shape(input)[1]  
    
    # get coordinates of array
    x = (fovx/2)*np.linspace(-1,1,x_size)
    y = (fovy/2)*np.linspace(-1,1,y_size)
    x_mesh, y_mesh = np.meshgrid(x,y)

    # plot projection line from pca eigenvector. Because of nearest neighbor interpolation, the axis 
    # vector end point is <size>-1. The number of samples is hard coded.
    if np.abs(xv) < np.abs(yv):
        y_line = np.linspace(0,y_size-1,nsample)
        x_line = xv / yv * y_line + y_size/2
        x_line = x_line + x_size/2 - x_line[int(nsample/2)]
    else:
        x_line = np.linspace(0,x_size-1,nsample)
        y_line = yv / xv * x_line + x_size/2
        y_line = y_line + y_size/2 - y_line[int(nsample/2)]
    
    xx_line = x_mesh[np.round(y_line).astype(int),np.round(x_line).astype(int)]
    yy_line = y_mesh[np.round(y_line).astype(int),np.round(x_line).astype(int)]
    
    d = np.sqrt(xx_line**2 + yy_line**2)
    
    # get minimum
    d_min = np.argwhere(np.min(d) == d)
    
    if np.size(d_min) > 1:
        d_min = np.asscalar(d_min[1])
    else:
        d_min = np.asscalar(d_min)
    
    d[:d_min] = -d[:d_min]
        
    # get autocorrelation
    array_acorr = get_acorr(input)
    acorr_line = array_acorr[np.round(y_line).astype(int),np.round(x_line).astype(int)]
    
    # fwhm
    acorr_line_max = np.argwhere(np.max(acorr_line) == acorr_line)
    if np.size(acorr_line_max) > 1:
        acorr_line_max = np.asscalar(acorr_line_max[0])
    else:
        acorr_line_max = np.asscalar(acorr_line_max)
    
    acorr_line_middle = acorr_line.copy()
    acorr_line_middle[acorr_line_middle > 0.5] = 0
    acorr_line_middle = np.argwhere(np.max(acorr_line_middle) == acorr_line_middle)
    
    if np.size(acorr_line_middle) > 1:
        acorr_line_middle = np.asscalar(acorr_line_middle[0])
    else:
        acorr_line_middle = np.asscalar(acorr_line_middle)

    fwhm_central = 2*np.abs(d[acorr_line_max] - d[acorr_line_middle])
    
    # spacing to neighbor
    peak = find_peaks(acorr_line, prominence=(p_min, p_max))[0]
    
    if len(peak) <= 1:
        d_neighbor = np.nan
        P_neighbor = np.nan
    else:
        peak_temp = peak - acorr_line_max
        peak_temp = np.abs(peak_temp)
        peak_temp = peak_temp.astype(float)
        peak_temp[peak_temp == 0] = np.NaN
        peak_temp = np.asscalar(np.argwhere(np.nanmin(peak_temp) == peak_temp)[0])
    
        d_neighbor_temp = peak[peak_temp]
        P_neighbor = acorr_line[d_neighbor_temp]
        d_neighbor = np.abs(d[acorr_line_max] - d[d_neighbor_temp])
    
    return fwhm_central, d_neighbor, P_neighbor, d, acorr_line