# Numerical methods course, AUT
# website: www.cemf.ir
'''
Newton's method to find the solution of a set of non-linear functions

inputs:
   fx: functions in the set
   x0: a vector of initial guess
   tol: relative tolerance 
   alpha: under-relaxation
output:
   [1]: solution of the set
   [2]: functions values at root
   [3]: number of iterations
   [4]: approximate error 
'''

from math import * 
from numpy import * #you should have numpy installed on your system

def newtonMethodSet2(fx, x0, tol = 1.0e-5, alpha = 1):

    #interval for numerical derivatives
    dx = max(0.00001, tol);
    nMax = 100; #maximum number of iterations

    # start newton loop
    for iter in range(nMax):
               
        (jac , f) = jacobianMat(fx, x0, dx);
        
        # check if x0 is the solution of the set
        ea = dot(f,f)
        if ea < tol:
            fx1 = f
            x1 = x0
            return (x0, f, iter, ea);
        
        delta = linalg.solve(jac, -f);
        x1 = x0 + alpha*delta
        
        ea = dot( delta, delta )/ max( max(x1), 1.0 )
        if ea  < tol:
            return (x1, fx(x1), iter, ea);
        
        x0 = x1;
    
    #no solution 
    raise Exception('newtonMethodSet2 did not converged in {} iteration'.format(i))

# Jacobian matrix for newton method (numerical derivatives)
def jacobianMat(fx, x, dx):

    n = len(x) #number of equations
    jac = zeros((n,n)) 
    fx0 = fx(x)
    
    for i in range(n):
        holdXi = x[i]
        x[i] = x[i]+ dx
        fx1 = fx(x)
        x[i] = holdXi
        jac[:,i] = (fx1-fx0)/dx #derivative with respect to xi

    return (jac, fx0);
