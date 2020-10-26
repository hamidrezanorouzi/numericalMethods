# Numerical methods course, AUT
# website: www.cemf.ir
def newton2(fx, x0 , tol = 1.0e-5):
    """
    Newton method to find one root of fx with the initial guess x0
    without the need for derivatives of fx
    DEFs:
    inputs:
        fx: function
        x0: initial guess
        tol: relative approximate tolerance 
    output:
        [1]: root
        [2]: function value at root
        [3]: number of iterations
        [4]: approximate error
    """


    nMax = 100 #maximum number of iterations
    dx = max(0.00001, tol) # interval for numerical drivatives

    fx0 = fx(x0)
    
    for i in range (1,nMax):

        fx0dx = fx(x0+dx)
        dfx0 = (fx0dx - fx0)/dx #forward derivatives
    
        x1 = x0 - fx0/dfx0
        fx1 = fx(x1);

        if abs(fx1) < tol :
            return (x1, fx1, i, abs(fx1));
        
        if abs(x1)>1.0e-15 :
            ea = abs((x1-x0)/x1)
        else:
            ea = abs(x1)

        if ea < tol :
            return (x1, fx1, i, ea);
        
        x0 = x1
        fx0 = fx1

    # no solution 
    raise Exception('newton2 did not converged in {} iteration'.format(nMax))

