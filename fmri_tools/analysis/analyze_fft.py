# -*- coding: utf-8 -*-

# external inputs
import numpy as np
from scipy.signal import find_peaks

# local inputs
from fmri_tools.utils import get_fft
    
    
def analyze_fft(input, fovx, fovy, xv, yv, f_cut=0.05, p_min=None, p_max=None, 
                nsample=1000):
    """
    This function computes the peak frequency and its corresponding power from a 
    one-sided power spectrum sampled on a projection lines of a 2d array using 
    nearest neighbor interpolation.
    Inputs:
        *input: input array (2D).
        *fovx: field of view in x-direction (mm).
        *fovy: field of view in y-firection (mm).
        *xv: x-coordinate of pca eigenvector.
        *yv: y-coordinate of pca eigenvector.
        *f_cut: cut off central spatial frequencies for peak detection.
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
    Last modified: 12-10-2020
    """

    # add one if nsample is an odd integer
    if np.mod(nsample,2) > 0:
        nsample += 1

    # get size of array
    x_size = np.shape(input)[0]
    y_size = np.shape(input)[1]
    
    # get coordinate of array
    kx = np.linspace(-1,1,x_size)*(x_size/(2*fovx))
    ky = np.linspace(-1,1,y_size)*(y_size/(2*fovy))
    ky_mesh, kx_mesh = np.meshgrid(ky,kx)

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
    
    kxx_line = kx_mesh[np.round(x_line).astype(int),np.round(y_line).astype(int)]
    kyy_line = ky_mesh[np.round(x_line).astype(int),np.round(y_line).astype(int)]
   
    # get final k-axis
    k_line = np.sqrt(kxx_line**2 + kyy_line**2)

    # get fourier spectrum        
    array_fft = get_fft(input)
    fft_line = array_fft[np.round(x_line).astype(int),np.round(y_line).astype(int)]
    
    # get one-sided spectrum
    arg_null = int(np.argwhere(k_line==np.min(k_line))[-1])
    fft_line = fft_line[arg_null:]
    k_line = k_line[arg_null:]
       
    # normalize by central frequency
    fft_line = fft_line / fft_line[0] * 100
    
    # first central k-space line for peak detection
    fft_cut = fft_line.copy()
    k_cut = k_line.copy()
    
    fft_cut = fft_cut[k_cut > f_cut]
    k_cut = k_cut[k_cut > f_cut]
    
    # find peaks
    peak = find_peaks(fft_cut, prominence=(p_min,p_max))[0]
    
    if len(peak) < 1:
        P_max = np.nan
        k_max = np.nan
    else:        
        P_max = fft_cut[peak[np.argwhere(fft_cut[peak] == np.max(fft_cut[peak]))[0]]]
        k_max = k_cut[peak[np.argwhere(fft_cut[peak] == np.max(fft_cut[peak]))[0]]]
    
    return k_max, P_max, k_line, fft_line
