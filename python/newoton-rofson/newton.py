import numpy as np


def newton(f, t0, initialy, h, num_steps):
    """
    Solve an ODE using the newton.

    Args:
        f: Function defining the ODE dy/dt = f(t, y).
        t0: Initial value of the independent variable.
        y0: Initial value of the dependent variable.
        h: Step size.
        num_steps: Number of steps to take.

    Returns:
        Tuple of arrays (t, y), where t contains the values of the independent variable
        and y contains the corresponding values of the dependent variable.
    """
    t = np.zeros(num_steps + 1)
    y = np.zeros(num_steps + 1)
    t[0] = t0
    y[0] = initialy

    for i in range(num_steps):
        k1 = h * f(t[i], y[i])
        k2 = h * f(t[i] + h / 2, y[i] + k1 / 2)
        k3 = h * f(t[i] + h / 2, y[i] + k2 / 2)
        k4 = h * f(t[i] + h, y[i] + k3)

        t[i + 1] = t[i] + h
        y[i + 1] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6

    return t, y
def my_ode(t, y):
    return t + y

t0 = 0.0  # Initial value of t
initialy = 1.0  # Initial value of y
h = 0.1   # Step size
num_steps = 10  # Number of steps

t, y = newton(my_ode, t0, initialy, h, num_steps)

# Print the results
for i in range(num_steps + 1):
    print(f"t = {t[i]}, y = {y[i]}")