def Heun(dydx, xspan, y0, n):

    """
    Modified Euler method to solve an ivp ODE

    parameter
    -----------------
        dydx: function to be integrated ATENTION!!! input lambda x,y: f even f had not y 
        xspan: [xi, xf] where xi and xf = initial and final values of independent variable
        y0: initial value of dependent variable
        n: number of steps

    retrun
    -----------------
        [1]: List of independent variable(t)
        [2]: List of solution for dependent variable(y)      
    """

    xi = xspan[0]; xf = xspan[1]
    if (xi>xf): raise exception("upper limit must be greater than lower")

    h = (xf - xi) / n

    #initial condition  
    y = [y0]; x = [xi]

    for i in range(n): #implement Heun's method
        
        iter = 100
        yf = GuessAndError (x[i], y[i], h, dydx, iter) # in The Heun method we take average
                                                            # between start and end slope of span
        x.append(x[i] + h)                                  # to have better estimate of y
        y.append(yf)                                        # y_end = y_start + h * (slope_start + slope_end)/2
                                                            # But when our dydx function is function of
                                                            # y and x we need to know y in end of span
                                                            # for calculate slope_end and we don't know it
                                                            # we can make an estimate with Euler (y_end = y_start + h * slope_start)
                                                            # But we can improve our estimate by using GuessAndError 
                                                            # because y_end = y_start + h * (slope_start + slope_end)/2
                                                            # and as we say slope_end is function of y_end and for that
                                                            # we can use GuessAndError to find better y_end

    return x, y        


def GuessAndError (xi, yi, h, dydx, iter):

    #Here xi and yi is not our initial x and y in whole span 
    slopei = dydx(xi, yi)

    yf0 = yi + slopei * h #first predict of yf
    
    tol = 1.0e-4;
    
    for i in range(iter):

        slopef = dydx(xi + h , yf0)
        yf = yi + (slopei+slopef)/2 * h

        Error = abs((yf-yf0)/ max(yf, tol)) #use max(ynew, tol) to avoids divided by zero

        if Error <= tol :
            print("shomare {}".format(i))
            return yf

        yf0 = yf


    raise Expection("solution of predictor-corrector did not converged in {} iterations".format(iter))
    return NaN
