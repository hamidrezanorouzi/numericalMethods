import numpy as np

def FalsePosition(f, xl, xu,  ea = 1e-3, Max_iter = 100):

    """
    False positon method for finding root in [xl, xu]
    f(xl) and f(xu) must have opposite sing
    ** This function wrote to find only one root in [xl, xu] **

    Note
    -----------
    in some cases (xr_new - xr_old) / xr_new < ea
    is True but f(xr_new) != 0 for That situation 
    function return Exception : False Position didn't converge

    parameters
    -----------

    f  : function that want to find it's root
    xl : lower limit
    xu : uper limit
    ea : percent of error default is 1e-3%
    Max_iter : maximume iteration for halving interval
    
    returns
    --------
    
    xr : the estimate of root
    
    """
    fxl = f(xl); fxu = f(xu)
    xr_old = xl

    if fxl * fxu > 0 : raise Exception("f(xl) and f(xu) must have opposite sings")

    for i in range(Max_iter):

        xr_new = xu - fxu * (xl - xu) / (fxl - fxu)
        fxr = f(xr_new)

        if fxl * fxr < 0 :
            xu = xr_new; fxu = fxr

        elif fxr * fxu < 0 :
            xl = xr_new; fxl = fxr

        elif abs(fxr) < ea:
            break

        Ea = abs ((xr_new - xr_old) / max(xr_new, 1.0e-15)) * 100

        if Ea < ea:
            break

        xr_old = xr_new

    if abs(fxr) > ea:
        raise Exception("False Position didn't converge")


    return xr_new


def ModifiedFalsePosition(f, xl, xu,  ea = 1e-3, Max_iter = 100):

    """ 
    for functions with significant curvature
    false position lead to poor convergence
    for that we have ModifiedFalsePosition
    that count how many time upper or lower bond
    didn't change and if any upper or lower bond
    stay fixed bigger that 2 time the function value
    in that bond is divided by 2 and again calculation is done

    parameters
    -----------

    f  : function that want to find it's root
    xl : lower limit
    xu : uper limit
    ea : percent of error default is 1e-3%
    Max_iter : maximume iteration for halving interval
    
    returns
    --------
    
    xr : the estimate of root

    """
    fxl = f(xl); fxu = f(xu)
    xr_old = xl

    l_iter = 0; u_iter = 0

    if fxl * fxu > 0 : raise Exception("f(xl) and f(xu) must have opposite sings")

    for i in range(Max_iter):

        xr_new = xu - fxu * (xl - xu) / (fxl - fxu)
        fxr = f(xr_new)

        if fxl * fxr < 0 :
            xu = xr_new; fxu = fxr

            u_iter = 0; l_iter = l_iter + 1
            if l_iter >= 2: fxl = fxl / 2

        elif fxr * fxu < 0 :
            xl = xr_new; fxl = fxr

            l_iter = 0; u_iter = u_iter + 1
            if u_iter >= 2: fxu = fxu / 2

        elif abs(fxr) < ea : break

        Ea = abs ((xr_new - xr_old) / max(xr_new, 1.0e-15)) * 100

        if Ea < ea:
            break

        xr_old = xr_new

    if abs(fxr) > ea:
        raise Exception("Modified False Position didn't converge")

    return xr_new