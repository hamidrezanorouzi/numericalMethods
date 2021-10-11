# Numerical methods course, AUT
# website: www.cemf.ir
import numpy as np
def triDiagonal(a, b, c, d):
    """
    Tomas method to solve a tri-diagonal set of linear equations 
    DEFs:
    inputs:
        a: lower diagonal elements [0:n]
        b: diagonal elements [0:n]
        c: upper diagonal elements [0:n]
        d: known coefficients [10n]
    output:
        X: solution of the set
    """
    #number of equations in the set
    n = len(b);
    
    cc = np.zeros((n,1));
    dd = np.zeros((n,1));
    
    #forward sweeping
    cc[0] = c[0]/b[0];
    for i in range(1,n-1):
        cc[i] = c[i]/(b[i]-a[i]*cc[i-1]);
    
    dd[0] = d[0]/b[0];    
    for i in range(1,n):
        dd[i] = (d[i]-a[i]*dd[i-1])/(b[i]-a[i]*cc[i-1]);
        
    # back substitution
    X = np.zeros((n,1));
    X[n-1] = dd[n-1];
    for i in range(n-2,-1,-1):
        X[i] = dd[i] - cc[i]*X[i+1];
    
    return X;    

