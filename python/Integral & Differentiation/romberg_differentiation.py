from differentiation import *
import numpy as np

def RombergDifferentiation(f, x, kmax = 4, err = 1.48e-8):

    """

    use RichardsonExtrapolation for calculatin differentiation

    parameters
    -------------
    f    : function we want it's derivative
    x    : where we want derivative to calculate
    kmax : order of error O(h**(kmax-1))
    err  : error

    return
    ------------

    diff, diff_matrix  : estimate of diff, the matrix of all order of error calculate
    
    """
    
    diff = np.zeros((kmax,kmax),float)

    for i in range(kmax):

        diff[i,0] = CenteralDifferentiation(f, x, h = 1/2**i, n = 1, error_order = 4)

        for k in range(i):

            diff[i, k+1] = RichardsonExtrapolation(diff[i][k], diff[i-1][k], k + 2)

    return diff[-1][-1], diff

def RichardsonExtrapolation(diff2, diff1, k):

    return (4**(k-1) * diff2 - diff1) / (4**(k-1) - 1)
