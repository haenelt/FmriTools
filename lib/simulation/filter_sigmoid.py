def filter_sigmoid(x, alpha=4):
    """
    Nonlinear sigmoid filter to adjust sharpness of ocular dominance. If sharpness factor is set to 
    zero, nothing is changed on the input image. The sigmoidal filter outputs an array with 
    intensity range [-1,1].
    Inputs:
        *x: input image.
        *alpha: sharpness factor.
            
    created by Daniel Haenelt
    Date created: 29-11-2018             
    Last modified: 09-01-2019
    """
    import numpy as np

    if alpha is not 0:            
        #return 1 / (1 + np.exp(-alpha*x))
        return 2 / ( 1 + np.exp(-alpha*x) ) -1
    else:
        return x