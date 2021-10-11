import numpy as np


def RK_Fehlberg(dydx, xspan, y0, n):

	"""
	in all RK2...4 methods we use constant step size
	for significant number of problem it can be limitation
	Runge-Kutta Fehlberg is part of ADAPTIVE RUNGE-KUTTA METHODS
	methods that decrease step size for the region of abrupt change
	and increase step size when changes is gradualy

	parameters
    ------------
    dydx  : function that want to find it's integral
    xspan : [a, b]
    y0    : initial condition
    n     : initial number of segment

    return
    -----------

    [1] : value of independent variable
    [2] : value of dependent variable
    
	"""

    [xi, xf] = xspan
    if (xi>xf): raise Exception("upper limit must be greater than lower")

    h = (xf - xi) / n

    #initial condition
    y = [y0]; x = [xi]

    x_now = xi


    Ea = []
    while xi<= x_now <= xf:

        y0_now = y[-1] # our initial conditon or (BC) for calcualte (yi+1) is y in previous step(last y)
        K1 = dydx(x_now, y0_now)
        K2 = dydx(x_now + h/5, y0_now + K1*h/5)
        K3 = dydx(x_now + 3*h/10, y0_now + 3/40 *K1*h + 9/40 * K2*h)
        K4 = dydx(x_now + 3/5*h, y0_now + 3/10 *K1*h - 9/10* K2*h + 6/5 * K3*h)
        K5 = dydx(x_now + h, y0_now - 11/54 * K1*h + 5/2 * K2*h - 70/27 * K3*h + 35/27 * K4*h)
        K6 = dydx(x_now + 7/8*h, y0_now + 1631/55296 * K1*h + 175/512 * K2*h + 575/13824 * K3*h
                 + 44275/110592 * K4*h + 253/4096 * K5*h)

        y_next_fourth = y0_now + (37/378*K1 + 250/621*K3 + 125/594*K4 + 512/1771*K6) * h
        y_next_fifth = y0_now + (2825/27648*K1 + 18575/48384*K3 + 13525/55296*K4 + 277/14336*K5 + K6/4) * h

        
        Ea_this_step = y_next_fifth - y_next_fourth 
        
        err_up = 1e-5
        err_down = 1e-8
        if abs(Ea_this_step) > err_up: # Time for decrease step size and back to previous step to recalculate
            h = h * abs(err_up/Ea_this_step)**0.25
            continue # recalculate with new h

        elif abs(Ea_this_step) < err_down: # Time for increase step size (for decrease computation) and back to previous step to recalculate
            h = h * abs(err_down/Ea_this_step)**0.2
            continue
        
        x_now = x_now + h
        x.append(x_now)
        y.append(y_next_fifth)

    return  x ,y
