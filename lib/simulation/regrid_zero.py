def regrid_zero_2d(data_array,Nx_new,Ny_new):
    """
    This function loads a two-dimensional numpy array and performs spatial interpolation using 
    zero padding in spatial frequency space. If the input array has an odd array size, the padding 
    size on the right size is decreased by one to match the output array size and leave the k-space
    center line at the same place as for even array sizes.
    Inputs:
        *data_array: input array.
        *Nx_new: resolution of output array in x-direction (only put even array size).
        *Ny_new: resolution of output array in y-direction (only put even array size).
    Outputs:
        *data_array_new: spatial interpolation of input array.
        
    created by Daniel Haenelt
    Date created: 07-01-2019      
    Last modified: 09-01-2019
    """
    import numpy as np
    from numpy.fft import fft2, ifft2, fftshift

    # size of input array
    data_size = np.shape(data_array)
    
    Nx_old = data_size[0]
    Ny_old = data_size[1]

    # zero padding size
    padx_left = int(Nx_new/2) - int(Nx_old/2)
    pady_left = int(Ny_new/2) - int(Ny_old/2)
    padx_right = int(Nx_new/2) - int(Nx_old/2)
    pady_right = int(Ny_new/2) - int(Ny_old/2)    
    
    # check if old data arrays is even/odd
    if Nx_old % 2:
        padx_right -= 1
   
    if Ny_old % 2:
        pady_right -= 1
    
    # Fourier transform of input data
    data_array_fft = fft2(data_array)
    data_array_fft = fftshift(data_array_fft)

    # pad zeros in spatial frequency space
    data_array_fft = np.pad(data_array_fft,
                            ((padx_left,padx_right),(pady_left,pady_right)),
                            'constant')    

    data_array_fft = fftshift(data_array_fft)

    # inverse Fourier transform of zero padded spatial frequency representation
    data_array_new = ifft2(data_array_fft)
    data_array_new = np.real(data_array_new)

    return data_array_new

def regrid_zero_1d(data_array,N_new):
    """
    This function loads a one-dimensional numpy array and performs spatial interpolation using 
    zero padding spatial frequency space. If the input array has an odd array size, the padding size
    on the right size is decreased by one to match the output array size and leave the k-space
    center line at the same place as for even array sizes.
    Inputs:
        *data_array: input array.
        *N_new: resolution of output array (only put even array size).
    Outputs:
        *data_array_new: spatial interpolation of input array.
        
    created by Daniel Haenelt
    Date created: 08-01-2019      
    Last modified: 09-01-2019
    """
    import numpy as np
    from numpy.fft import fft, ifft, fftshift

    # size of input array
    data_size = np.shape(data_array)
    
    N_old = data_size[0]

    # zero padding size
    pad_left = int(N_new/2) - int(N_old/2)
    pad_right = int(N_new/2) - int(N_old/2)
    
    # check if old or new data arrays are even/odd
    if N_old % 2:
        pad_right -= 1
    
    # Fourier transform of input data
    data_array_fft = fft(data_array)
    data_array_fft = fftshift(data_array_fft)

    # pad zeros in spatial frequency space
    data_array_fft = np.pad(data_array_fft, (pad_left,pad_right), 'constant')    
            
    data_array_fft = fftshift(data_array_fft)

    # inverse Fourier transform of zero padded spatial frequency representation
    data_array_new = ifft(data_array_fft)
    data_array_new = np.real(data_array_new)

    return data_array_new