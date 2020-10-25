# Numerical methods course, AUT
# website: www.cemf.ir
def secant(fx, xs, xe, tol = 1.0e-5):
    """
    Secant method to find one root of fx in [xs,xe]
    DEFs:
    inputs:
        fx: function
        xs: start point of the search interval
        xe: end point of the search interval
    output:
        [1]: root
        [2]: function value at root
        [3]: number of iterations
        [4]: approximate error
    """

    fs = fx(xs)
    fe = fx(xe)
    xz_old = xs;


    for i in range(1,1000):

        xz = xe - (xe-xs)*fe/(fe-fs)
        fz = fx(xz)

        if abs(fz) < tol :
            return (xz, fz, i, fz);

        if fz*fe <0 :
            xs = xz
            fs = fz
        else:
            xe = xz
            fe = fz

        if abs(xz) > 1.0e-15:
            ea = abs((xz-xz_old)/xz)
        else:
            ea = abs(xz)

        if ea < tol:
            return (xz, fz, i, ea);

        xz_old = xz

    raise Exception('secant did not converged in {} iteration'.format(i))

        
            
