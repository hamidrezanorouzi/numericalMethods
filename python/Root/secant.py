def Secant(fx, x0, x1, tol = 1.0e-5, Max_iter = 100):

    """
    Secant method different from newton raphson
    in calculating derivative of fx in newton raphson
    we have fx' but here we use 2 inital guess 
    and approximate derivative by a straight line 
    between two points(backward finite divided difference)

    parameters
    ------------
    fx  : function that want to find it's root
    x0  : zero inital guess
    x1  : first inital guess
    tol : error default is 1e-5
    Max_iter : maximume iteration for loop
    
    returns
    --------
    x2 : root guess
    fx2: function value in x2
    i  : number of iterations

    """

    fx0 = fx(x0)
    fx1 = fx(x1)

    for i in range (Max_iter):

        x2 = x1 - (fx1* (x0 - x1)) / (fx0 - fx1)

        fx2 = fx(x2)

        if abs(fx2) < tol:

            return (x2, fx2, i);

        x0 = x1; x1 = x2
        fx0 = fx1; fx1 = fx2

    
    raise Exception('Secant did not converged in {} iteration'.format(nMax))


def ModifiedSecant(fx, x0, tol = 1.0e-5, Max_iter = 100):

    
    """
    Modified Secant method instead of input 
    2 guess get 1 guess and approximate derivative 
    by tiny step size max(0.00001, tol)

    parameters
    ------------
    fx  : function that want to find it's root
    x0  : inital guess
    tol : error default is 1e-5
    Max_iter : maximume iteration for loop
    
    returns
    --------
    x1 : root guess
    fx1: function value in x2
    i  : number of iterations

    Noted 
    ---------
    If dx is too small, the method can be swamped by round-off error 
    caused by subtractive cancellation in the denominator
    If it is too big, the technique can become inefficient and even divergent.
    However, if chosen correctly, it provides a nice alternative for cases where evaluating
    the derivative is diffi cult and developing two initial guesses is inconvenient.

    """
    
    dx = max(0.00001, tol) # interval for numerical drivatives

    fx0 = fx(x0)
    
    for i in range (Max_iter):


        x1 = x0 - (dx * fx(x0)) / (fx(x0+dx) - fx(x0))
        
        fx1 = fx(x1);

        if abs(fx1) < tol :

            return (x1, fx1, i);
        
        x0 = x1
        fx0 = fx1
    
    raise Exception('Modified Secant did not converged in {} iteration'.format(Max_iter))
