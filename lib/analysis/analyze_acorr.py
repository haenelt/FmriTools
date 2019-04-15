def analyze_acorr(input, fovx, fovy, xv, yv, p_min=0.01, p_max=None, nsample=1000):
    """
    This function computes the normalized autocorrelation (NAC) from a 2D input array and estimates 
    the width of the central peak and the distance to its first neighbor peak along a defined 
    projection lines sampled with nearest neighbor interpolation. The width of the central peak is 
    defined as the width at zero point.
    Inputs:
        *input: input array (2D).
        *fovx: field of view in x-direction (mm).
        *fovy: field of view in y-firection (mm).
        *xv: x-coordinate of pca eigenvector.
        *yv: y-coordinate of pca eigenvector.
        *p_min: minimum prominence for peak detection.
        *p_max: maximum prominence for peak detection.
        *nsample: number of sampling points along projection line.
    Outputs:
        *fwhm_central: fwhm of nac central peak in mm along projection axis.
        *P_neighbor: power of first neighbor peak along projection axis.
        *d_neighbor: distance to first neighbor peak in mm along projection axis.
        *d: lag in mm along projection axis.
        *acorr_line: nac along projection axis.
        
    created by Daniel Haenelt
    Date created: 14-04-2019
    Last modified: 15-04-2019
    """
    import numpy as np
    import copy
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

    # get projection line coordinates from pca eigenvector. Because we interpolate using nearest
    # neighbors the axis enc point is <size>-1. 
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
    
    # get final distances
    d = np.sqrt(xx_line**2 + yy_line**2)
    
    # get minimum to shift mid point to origin
    d_min = np.argwhere(np.min(d) == d)
    
    if np.size(d_min) > 1:
        d_min = int(d_min[-1])
    else:
        d_min = int(d_min)
    
    # all coordinates left from origin get negative signs
    d[:d_min] = -d[:d_min]
        
    # get autocorrelation
    array_acorr = get_acorr(input)
    acorr_line = array_acorr[np.round(y_line).astype(int),np.round(x_line).astype(int)]
    
    # fwhm
    acorr_line_max = np.argwhere(np.max(acorr_line) == acorr_line)
    if np.size(acorr_line_max) > 1:
        acorr_line_max = int(acorr_line_max[0])
    else:
        acorr_line_max = int(acorr_line_max)
    
    acorr_line_middle = copy.deepcopy(acorr_line_max)
    while True:
        acorr_line_middle += 1
        if acorr_line_middle > nsample-1:
            acorr_line_middle = np.nan
            break
        elif acorr_line[acorr_line_middle] < 0:
            acorr_line_middle -= 1
            break
    
    # compute fwhm if middle point is found
    if ~np.isnan(acorr_line_middle):   
        fwhm_central = 2*np.abs(d[acorr_line_max] - d[acorr_line_middle])
    else:
        fwhm_central = np.nan
        
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
        peak_temp = int(np.argwhere(np.nanmin(peak_temp) == peak_temp)[0])
    
        P_neighbor = acorr_line[peak[peak_temp]]
        d_neighbor = np.abs(d[acorr_line_max] - d[peak[peak_temp]])
    
    return fwhm_central, d_neighbor, P_neighbor, d, acorr_line