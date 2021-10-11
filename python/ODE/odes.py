import numpy as np

def ODEs(List_of_dydx, xspan, List_of_y0, n):

    """
    solve system of ODEs with fourth runge katta

    parameter
    -----------------
        List_of_dydx: all functions to be integrated simultaneously( lambda x,y1,y2,...yn : [f1, f2, ..., fn] )
        xspan: [xi, xf] where xi and xf = initial and final values of independent variable
        y0: initial values of dependent variable
        n: number of steps

    return
    -----------------
        [1]: List of independent variable (x)
        [2]: List of solution for dependent variables (y) (n+1)*2 matrix in row we value 
             for dependent variables in x 
    """

    [a, b] = xspan
    if (a>b): raise Exception("upper limit must be greater than lower")

    h = (b - a) / n

    #initial condition
    y = [List_of_y0] ; x = [a]

    for this_step in range(n):

        y_this_step = []
        span = [x[this_step], x[this_step]+h]
        
        x_in_this_step, y_in_this_step = RK4_systems(List_of_dydx, span, y[this_step], 1) 

        x.append(x_in_this_step[1])
        y.append(y_in_this_step[1])

    return x, y


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

    #initial condition
    y = [y0]; x = [xi]

    for i in range(n):
        
        x.append(x[i] + h)
        yi = np.array(y[i])

        K1 = np.array(List_of_dydx(x[i], *yi))
        K2 = np.array(List_of_dydx(x[i] + h/2, *yi + K1 * h/2))
        K3 = np.array(List_of_dydx(x[i] + h/2, *yi + K2 * h/2))
        K4 = np.array(List_of_dydx(x[i] + h, *yi + K3 * h))

        y.append((yi + h* (K1 + 2*K2 + 2*K3 + K4)/6).tolist())

    return x, y

