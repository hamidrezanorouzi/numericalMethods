import numpy as np

def Trapezoidal(f, span, n = 1):

    """
    Trapezoidal for evaluating
    itegral of function it'sn't for 
    tabulated data

    parameters
    ------------
    f    : function in integral
    span : [a, b]
    n    : number of parts

    return 
    ------------
    I : estimated value of integral

    """
    [a, b] = span
    if ( a > b ): raise Exception("upper limit must be greater than lower")

    x = np.linspace(a, b, n+1)
    fx = f(x)

    I = (b - a) * (fx[0] + 2*sum(fx[1:-1]) + fx[-1]) / (2*n)
    return I
