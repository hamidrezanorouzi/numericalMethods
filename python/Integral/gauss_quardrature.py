import numpy as np

weighting_factor = [[1,1],
                    [0.5555556, 0.8888889, 0.5555556],
                    [0.3478548, 0.6521452, 0.6521452, 0.3478548],
                    [0.2369269, 0.4786287, 0.5688889, 0.4786287, 0.2369269],
                    [0.1713245, 0.3607616, 0.4679139, 0.4679139, 0.3607616, 0.1713245]]

function_argument = [[-0.577350269, 0.577350269],
                     [-0.774596669, 0, 0.774596669],
                     [-0.861136312, -0.339981044, 0.339981044,0.861136312],
                     [-0.906179846, -0.538469310, 0 , 0.538469310, 0.906179846],
                     [-0.932469514, -0.661209386, -0.238619186, 0.238619186, 0.661209386, 0.932469514]]


def functionTransform(f, span):

    """
        transform the integration interval without changing 
        the value of the integral. I(f[a,b]) = I (fT[-1,1])

        f span = [a,b] ---> fTransfer span = [-1,1]

        x = (b+a + (b-a)* xd)/2
        dx = (b-a)/2 * dxd
    """
    [a, b] = span
    f_Trans = lambda xd : f(((b+a) + (b-a)*xd)/2) * (b-a)/2
    
    return f_Trans


def GaussLegendre(f, span, n = 2):

    """
    I = c0*fT(x0) + c1*fT(x1) + ... + cn-1 * fT(xn-1)
    values for c's and x's for up to and including the six 
    point formula

    parameter
    -----------------

    f    : function in integral
    span : span of integral
    n    : number of point in GaussLegendre formula

    return
    -----------------

    I    : integral estimation

    """

    fT = functionTransform(f, span)
    row_of_metrix = n-2

    c = np.array(weighting_factor[row_of_metrix])
    x = np.array(function_argument[row_of_metrix])

    I = sum(c * fT(x))

    return I