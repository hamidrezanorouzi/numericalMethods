import numpy as np

def NewtonRaphson (fx, dfx, x0, tol = 1e-5):

    """
    relative error isn't write

    """

    Max_iter = 100;
    fx0 = fx(x0);

    for i in range(Max_iter):

        x1 = x0 - fx0 / max(dfx(x0), tol); # use max(dfx(x0), tol) for preventing dividing by zero
        fx1 = fx(x1);
        
        if abs(fx1) <= tol:

            return (x1, fx1, i);
            
        x0 = x1;
        fx0 = fx1;

    raise Exception("Newton Raphson didn't converge in {} iterations".format(Max_iter))

