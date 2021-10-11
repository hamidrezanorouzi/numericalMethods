import numpy as np

def FixedPoint(g, x0, ea = 1e-3, Max_iter = 100):

    """
    x = g(x)
    x_i+1 = g(x_i)

    parameters
    -----------

    g  : function that want to find it's root
    x0 : initial guess
    ea : percent of error default is 1e-3%
    Max_iter : maximume iteration
    
    returns
    --------
    
    xr : the estimate of root

    """
    xr = x0

    for i in range(Max_iter):

        xr_old = xr
        xr = g(xr_old)
        
        Ea = abs((xr - xr_old) / max(xr, 1.0e-15))
        
        if Ea < ea : 

            return xr, i

    raise Exception("Fixed Point didn't converge to {} error in {} iteration".format(ea, i))
    
