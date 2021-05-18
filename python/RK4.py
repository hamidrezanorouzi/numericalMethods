# Numerical methods course, AUT
# website: www.cemf.ir

import numpy as np 

def RK4(dydt, tspan, y0, n):

    """
    fourth order rk method to solve an ivp ODE
    
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
    y = []
    t = []

    #initial condition
    y.append(y0); t.append(ti)

    for i in range(n):

        t.append(t[i] + h)
        yi = np.array(y[i]) # beacuase with List type we can't Done Math calculation

        K1 = np.array(dydt(t[i], *yi))                  # This guy(*) convert List to argument 
        K2 = np.array(dydt(t[i] + h/2, *(yi + K1*h/2))) # beacuse our dydt lambda function
        K3 = np.array(dydt(t[i] + h/2, *(yi + K2*h/2))) # need n+1 inputs.
        K4 = np.array(dydt(t[i] + h, *(yi + K3*h)))     
                                                        
        

        y.append( yi + (K1 + 2*K2 + 2*K3 + K4)*h/6)

    return t, y