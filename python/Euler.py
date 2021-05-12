# Numerical methods course, AUT
# website: www.cemf.ir

def Euler(dydt, tspan, y0, n):
    """
    Euler Method for solve an ivp ODEs
    DEFs:
    inputs:
        dydt: function to be integrated
        tspan: [ti, tf] where ti and tf = initial and final values of independent variable
        y0: initial value of dependent variable
        n: number of steps
    outputs:
        [1]: vector of independent variable (t)
        [2]: vector of solution for dependent variable (y)

    """
    ti = tspan[0]; tf = tspan[1]

    if(tf < ti):
        raise Exception("upper limit must be greater than lower")

    h = (tf-ti)/n
    y = []; t=[]

    #initial condition
    y.append(y0); t.append(ti)

    for i in range(n):#implement Euler's method

        t.append(t[i] + h)
        y.append( y[i] + h * dydt(t[i], y[i]))
    
    return t, y
