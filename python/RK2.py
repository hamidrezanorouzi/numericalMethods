# Numerical methods course, AUT
# website: www.cemf.ir

import numpy as np

def RK2(dydt, tspan, y0, n):

    """
    Second order rk method to solve an ivp ODE
    
    dy1/dt = f(t, y1, ..., yn)
    dy2/dt = f2(t, y1, ..., yn)
    .
    .
    .
    dyn/dt = fn(t, y1, ..., yn)

    DEFs:
    inputs:
        dydt: function to be integrated ( dydt = lambda t, y1,...,yn: [f, f2, ..., fn] )
        tspan: [ti, tf] where ti and tf = initial and final values of independent variable
        y0: initial values of dependent variable
        n: number of steps

    output:
        [1]: vector of independent variable(t)
        [2]: vector of solution for dependent variable(y)      
    """

    ti = tspan[0]; tf = tspan[1]

    if (ti>tf):
        raise exception("upper limit must be greater than lower")

    h = (tf - ti) / n
    y = []
    t = []

    #initial condition
    y.append(y0); t.append(ti)

    for i in range(n):#implement rk method
        
        t.append(t[i] + h)
        yi = np.array(y[i]) # beacuase with List type we can't Done Math calculation

        K1 = np.array(dydt(t[i], *y[i]))            # This guy(*) convert List to argument 
        K2 = np.array(dydt(t[i]+h, *(yi+ h*K1)))    # beacuse our dydt, lambda function
                                                    # need n+1 inputs.
        y.append(y[i] + 0.5*(K1+K2)*h)

    return t, y