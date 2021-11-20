import numpy as np

def RK4_systems(List_of_dydx, xspan, y0, n):

    """
    fourth order Runge-Kutta method to solve an ivp ODEs
    DEFs:

    parameter
    -----------------
        List_of_dydx: functions to be integrated ( List_of_dydx = lambda t, y1,...,yn: [f, f2, ..., fn] )
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
    numberofequation = len(y0)
    x = np.zeros(n+1) ; y = np.zeros((numberofequation, n+1))

    #initial condition
    y[:, 0] = y0; x[0] = xi

    for i in range(n):
        
        x[i+1] = x[i] + h

        K1 = np.array(List_of_dydx(x[i], *y[:,i]))
        K2 = np.array(List_of_dydx(x[i] + h/2, *y[:,i] + K1 * h/2))
        K3 = np.array(List_of_dydx(x[i] + h/2, *y[:,i] + K2 * h/2))
        K4 = np.array(List_of_dydx(x[i] + h, *y[:,i] + K3 * h))

        y[:, i+1] = y[:, i] + h* (K1 + 2*K2 + 2*K3 + K4)/6

    return x, y

