# Numerical methods course, AUT
# website: www.cemf.ir

def modifiedEulerMethod(dydt, tspan, y0, n):

    """
    Modified Euler method to solve an ivp ODE
    DEFs:
    inputs:
        dydt: function to be integrated ATENTION!!! input lambda t,y: f even f had not y 
        tspan: [ti, tf] where ti and tf = initial and final values of independent variable
        y0: initial value of dependent variable
        n: number of steps

    output:
        [1]: vector of independent variable(t)
        [2]: vector of solution for dependent variable(y)      
    """

    ti = tspan[0]; tf = tspan[1]

    if (ti>tf):
        raise exception("upper limit must be greater than lower")

    h = (tf - ti) / n
    tol = 1.0e-4;
    y = []
    t = []

    #initial condition  
    y.append(y0); t.append(ti)

    for i in range(n): #implement Heun's method

        t.append(t[i] + h)
        fxi = dydt(t[i], y[i])

        #predictor equation 
        yold = y[i] + h*fxi

        converge = False
        for iter in range(1, 1000):

            fxi1 = dydt(t[i+1], yold)
            #Corrector  equation
            ynew = y[i] + h * (fxi + fxi1)/2
            err = abs((ynew-yold)/ max(ynew, tol)) #use max(ynew, tol) to avoids divided by zero
            if (err < tol):

                converge = True
                y.append(ynew)
                break

            else:
                yold = ynew

        if converge == False :
            raise Expection("solution of predictor-corrector did not converged in 1000 iterations")

    return t, y