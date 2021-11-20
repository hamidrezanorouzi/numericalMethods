# Numerical methods course, AUT
# website: www.cemf.ir

import numpy as np

def RK2(dydx, xspan, y0, n):

    """
    Second order Runge-kutta method to solve an ivp ODE
    RK2 is Huen method without iteration(GuessAndError for find better y)
    
    parameter
    -----------------
        dydx: function to be integrated ( ATENTION!!! input lambda x,y: f even f had not y  )
        xspan: [xi, xf] where xi and xf = initial and final values of independent variable
        y0: initial values of dependent variable
        n: number of steps

    return
    -----------------
        [1]: vector of independent variable(x)
        [2]: vector of solution for dependent variable(y1,...,yn)      
    """

    xi = xspan[0]; xf = xspan[1]

    if (xi>xf): raise exception("upper limit must be greater than lower")

    h = (xf - xi) / n

    x = np.zeros(n+1) ; y = np.zeros(n+1)

    #initial condition
    y[0] = y0; x[0] = xi

    for i in range(n):
        
        K1 = dydx(x[i], y[i])     
        K2 = dydx(x[i]+h, y[i]+ h*K1)
                         
        x[i+1] = x[i] + h
        y[i+1] = y[i] + 0.5*(K1+K2)*h

    return x, y

def RK3(dydx, xspan, y0, n):

    """
    third order Runge-kutta method to solve an ivp ODE

    parameter
    -----------------
        dydx: function to be integrated ( ATENTION!!! input lambda x,y: f even f had not y  )
        xspan: [xi, xf] where xi and xf = initial and final values of independent variable
        y0: initial values of dependent variable
        n: number of steps

    return
    -----------------
        [1]: vector of independent variable(x)
        [2]: vector of solution for dependent variable(y1,...,yn)      
    """

    xi = xspan[0]; xf = xspan[1]

    if (xi>xf): raise exception("upper limit must be greater than lower")

    h = (xf - xi) / n
    x = np.zeros(n+1) ; y = np.zeros(n+1)

    #initial condition
    y[0] = y0; x[0] = xi

    for i in range(n):
        
        
        
        K1 = dydx(x[i], y[i])     
        K2 = dydx(x[i]+ h/2, y[i]+ h/2*K1)
        K3 = dydx(x[i] + h, y[i] - K1*h + 2*K2*h)
        
        x[i+1] = x[i] + h                                           
        y[i+1] = y[i] + (K1 + 4*K2 + K3)/6 * h

    return x, y

def RK4(dydx, xspan, y0, n):

    """
    fourth order Runge-Kutta method to solve an ivp ODE

    parameter
    -----------------
        dydx: function to be integrated ( dydx = lambda t, y1,...,yn: [f, f2, ..., fn] )
        tspan: [ti, tf] where ti and tf = initial and final values of independent variable
        y0: initial value of dependent variable
        n: number of steps

    return
    -----------------
        [1]: vector of independent variable(t)
        [2]: vector of solution for dependent variable(y)      
    """
    xi = xspan[0]; xf = xspan[1]

    if (xi>xf): raise exception("upper limit must be greater than lower")

    h = (xf - xi) / n
    x = np.zeros(n+1) ; y = np.zeros(n+1)

    #initial condition
    y[0] = y0; x[0] = xi

    for i in range(n):
       
        K1 = dydx(x[i], y[i])     
        K2 = dydx(x[i] + h/2, y[i] + K1 * h/2)
        K3 = dydx(x[i] + h/2, y[i] + K2 * h/2)
        K4 = dydx(x[i] + h, y[i] + K3 * h)

        x[i+1] = x[i] + h
        y[i+1] = y[i] + h* (K1 + 2*K2 + 2*K3 + K4)/6

    return x, y


def RK5(dydx, xspan, y0, n):

    """
    fifth order Runge-kutta method(Butcher) to solve an ivp ODE

    parameter
    -----------------
        dydx: function to be integrated ( ATENTION!!! input lambda x,y: f even f had not y  )
        xspan: [xi, xf] where xi and xf = initial and final values of independent variable
        y0: initial values of dependent variable
        n: number of steps

    return
    -----------------
        [1]: vector of independent variable(x)
        [2]: vector of solution for dependent variable(y1,...,yn)      
    """

    xi = xspan[0]; xf = xspan[1]

    if (xi>xf): raise Exception("upper limit must be greater than lower")

    h = (xf - xi) / n
    x = np.zeros(n+1) ; y = np.zeros(n+1)

    #initial condition
    y[0] = y0; x[0] = xi

    for i in range(n):
        
        K1 = dydx(x[i], y[i])     
        K2 = dydx(x[i] + h/4, y[i] + K1*h/4)
        K3 = dydx(x[i]+ h/4, y[i]+ K1*h/8 + K2*h/8)
        K4 = dydx(x[i] + h/2, y[i] - K2*h/2 + K3*h)
        K5 = dydx(x[i] + 3/4*h, y[i] + 3/16*K1*h + 9/16*K4*h)
        K6 = dydx(x[i] + h, y[i] - 3/7*K1*h + 2/7*K2*h + 12/7*K3*h -12/7*K4*h + 8/7*K5*h)
        
        x[i+1] = x[i] + h
        y[i+1] = y[i] + (7*K1 + 32*K3 + 12*K4 + 32*K5 + 7*K6)/90 * h

    return x, y