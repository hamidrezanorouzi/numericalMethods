# Numerical methods course, AUT
# website: www.cemf.ir
def bisection(fx, a, b, tol = 1.0e-5):
    """
    bisection method to one find root of fx in [a,b]
    DEFs:
    inputs:
        fx: function
        a: start point of search interval
        b: end point of search interval
    output:
        [1]: root
        [2]: function value at root
        [3]: number of iterations
        [4]: approximate error
    """
    fa = fx(a)
    fb = fx(b)
    c_old = a;

    for i in range(1,1000):

        c = (a+b)/2
        fc = fx(c)

        if abs(fc) < tol :
            return (c, fc, i, fc);

        if fc*fb <0 :
            a = c
            fa = fc
        else:
            b = c
            fb = fc

        if abs(c) > 1.0e-15:
            ea = abs((c-c_old)/c)
        else:
            ea = abs(c)

        if ea < tol:
            return (c, fc, i, ea);

        c_old = c

    raise Exception('bisection did not converged in {} iteration'.format(i))

        
            
