# -*- coding: utf-8 -*-

# python standard library inputs
import copy

# external inputs
import numpy as np
from scipy.signal import find_peaks

# local inputs
from fmri_tools.utils.get_acorr import get_acorr
    
    
def analyze_acorr(input, fovx, fovy, xv, yv, p_min=0.01, p_max=0.5, 
                  nsample=1000):
    """ Analyze Autocorrelation
    
    This function computes the normalized autocorrelation (NAC) from a 2D input 
    array and estimates the width of the central peak and the distance to its 
    first neighbor peak along a defined projection lines sampled with nearest 
    neighbor interpolation. The width of the central peak is defined as the 
    width at zero point.

    Parameters
    ----------
    input : ndarray
        2D input array.
    fovx : float
        Field of view in x-direction (mm).
    fovy : float
        Field of view in y-firection (mm).
    xv : float
        x-coordinate of pca eigenvector.
    yv : float
        y-coordinate of pca eigenvector.
    p_min : float, optional
        Minimum prominence for peak detection. The default is 0.01.
    p_max : float, optional
        Maximum prominence for peak detection. The default is 0.5.
    nsample : int, optional
        Number of sampling points along projection line. The default is 1000.

    Returns
    -------
    fwhm_central : float
        FWHM of NAC central peak in mm along projection axis.
    d_neighbor : float
        Distance to first neighbor peak in mm along projection axis.
    P_neighbor : float
        Power of first neighbor peak along projection axis.
    d : float
        Lag in mm along projection axis.
    acorr_line : ndarray
        NAC along projection axis.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 14-04-2019
    Last modified: 13-10-2020
    
    """
    
    # add one if nsample is an odd integer
    if np.mod(nsample,2) > 0:
        nsample += 1

    # get size of array
    x_size = np.shape(input)[0]
    y_size = np.shape(input)[1]  
    
    # get coordinates of array
    x = (fovx/2)*np.linspace(-1,1,x_size)
    y = (fovy/2)*np.linspace(-1,1,y_size)
    y_mesh, x_mesh = np.meshgrid(y,x)

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
    
    # check that no line steps over array border
    if np.any(y_line > y_size-1):
        x_line[y_line > y_size-1] = np.nan
        y_line[y_line > y_size-1] = np.nan
        x_line = x_line[~np.isnan(x_line)]
        y_line = y_line[~np.isnan(y_line)]
    
    if np.any(y_line < 0):
        x_line[y_line < 0] = np.nan
        y_line[y_line < 0] = np.nan
        x_line = x_line[~np.isnan(x_line)]
        y_line = y_line[~np.isnan(y_line)]
    
    if np.any(x_line > x_size-1):
        y_line[x_line > x_size-1] = np.nan
        x_line[x_line > x_size-1] = np.nan
        x_line = x_line[~np.isnan(x_line)]
        y_line = y_line[~np.isnan(y_line)]
    
    if np.any(y_line < 0):
        y_line[x_line < 0] = np.nan
        x_line[x_line < 0] = np.nan
        x_line = x_line[~np.isnan(x_line)]
        y_line = y_line[~np.isnan(y_line)]
    
    xx_line = x_mesh[np.round(x_line).astype(int),np.round(y_line).astype(int)]
    yy_line = y_mesh[np.round(x_line).astype(int),np.round(y_line).astype(int)]
    
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
    acorr_line = array_acorr[np.round(x_line).astype(int),np.round(y_line).astype(int)]
        
    # fwhm
    acorr_line_max = np.argwhere(np.max(acorr_line) == acorr_line)
    if np.size(acorr_line_max) > 1:
        acorr_line_max = int(acorr_line_max[0])
    else:
        acorr_line_max = int(acorr_line_max)
    
    # FWHM will be defined at half maximum. N.B., this underestimates the columnar width in case of 
    # pure sinusoidal oscillation where the width would be determined by taking the FWHM at zero
    acorr_line_middle = copy.deepcopy(acorr_line_max)
    while True:
        acorr_line_middle += 1
        if acorr_line_middle > len(acorr_line)-1:
            acorr_line_middle = np.nan
            break
        elif acorr_line[acorr_line_middle] < 0.5:
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
