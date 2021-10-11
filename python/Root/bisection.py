def Bisection(f, xl, xu, ea = 1e-3, Max_iter = 100):

    """
    Bisection method for finding root in [xl, xu]
    f(xl) and f(xu) must have opposite sing
    ** This function wrote to find only one root in [xl, xu] **

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

    if fxl * fxu > 0 : raise Exception("f(xl) and f(xu) must have opposite sings")

    for i in range(Max_iter):

        xr = (xl + xu) / 2
        fxr = f(xr)

        if fxl * fxr < 0 :
            xu = xr; fxu = fxr

        elif fxr * fxu < 0 :
            xl = xr; fxl = fxr

        Ea = abs((xu - xl) / (xu + xl)) * 100

        if abs(fxr) < ea or Ea < ea:
            break

    return xr
