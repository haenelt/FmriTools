def get_white_2d(nx, ny, mu, sigma):
    """
    Creates a 2D array with Gaussian noise.
    Input:
        *nx: array size in x-direction.
        *ny: array size in y-direction.
        *mu: mean of the Gaussian distribution.
        *sigma: standard deviation of the Gaussian distribution.
    Output:
        *img: white noise array.

    created by Daniel Haenelt
    Date created: 04-01-2019
    Last modified: 07-01-2019
    """
    import numpy as np

    img = np.random.normal(mu,sigma,nx*ny)
    img = np.reshape(img,(nx,ny))
    
    return img


def get_white_1d(n, mu, sigma):
    """
    Creates a 1D array with Gaussian noise.
    Input:
        *n: array size.
        *mu: mean of the Gaussian distribution.
        *sigma: standard deviation of the Gaussian distribution.
    Output:
        *img: white noise array.

    created by Daniel Haenelt
    Date created: 08-01-2019
    Last modified: 08-01-2019
    """
    import numpy as np

    img = np.random.normal(mu,sigma,n)
    
    return img
