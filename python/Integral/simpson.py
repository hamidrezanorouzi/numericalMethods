import numpy as np

def simpson1_3(f, span, n = 2):

    """
    1/3 simpson rule

    Notes
    -------------
    n must be 2*N

    parameters
    -------------
    f    : function in integral
    span : [a, b] span of integral
    n    : numbers of section

    return
    -------------

    I : estimation of integral
    """
    [a, b] = span
    if ( a > b ): raise exception("upper limit must be greater than lower")
    if (n % 2 != 0) : raise exception("simpson1/3 is only for odd point or even segments")

    
    f0 = f(a); fn = f(b)

    x = np.linspace(a, b, n+1)

    I = ( b - a ) * (f0 + 4*sum(f(x[1:n:2])) + 2*sum(f(x[2:n-1:2])) + fn) / (3*n)
    return I


def simpson3_8(f, span , n = 3):

    """
    3/8 simpson rule

    Notes
    -------------
    n must be 3*N

    parameters
    -------------
    f    : function in integral
    span : [a, b] span of integral
    n    : numbers of section
    
    return
    -------------

    I : estimation of integral
    """

    [a, b] = span
    if ( a > b ): raise Exception("upper limit must be greater than lower")
    if (n % 3 != 0) : raise Exception("simpson3/8 is only for 3 multiple parts 3/6/9/...")

    
    f0 = f(a); fn = f(b)
    x = np.linspace(a, b, n+1)

    I = 3/8 * (b - a)/n * (f0 + 2*sum(f(x[3:n:3])) + 3*(sum(f(x[1:n:3])) + sum(f(x[2:n:3]))) + fn)
    return I
