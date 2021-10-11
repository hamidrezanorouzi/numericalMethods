from trapezoidal import *

def Romberg(f, span, kmax = 8):

    """
    parameters
    ------------
    f    : function in integral
    span : [a, b]
    kmax : size of square matrix
        
    return 
    ------------
    I : estimated value of integral with error order O(h^(2*k_max))
    I-matrix : the matrix of integrals

    """
    I = np.zeros((kmax,kmax),float)

    for i in range(kmax):

        I[i,0] = Trapezoidal(f, span, n = 2**i)

        for k in range(i):

            I[i, k+1] = RichardsonExtrapolation(I[i][k], I[i-1][k], k + 2)

    return I[-1][-1], I

def RichardsonExtrapolation(I2, I1, k):


    """
    Error-correction techniques are available to improve
    the results of numerical integration on the basis of the integral estimates themselves.
    Generally called Richardson’s extrapolation, these methods use two estimates of an integral
    to compute a third, more accurate approximation.

    """

    return (4**(k-1) * I2 - I1) / (4**(k-1) - 1)