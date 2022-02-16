import numpy as np

def LinearRegression(X, Y):

	"""

	linear regression of given data(X, Y) 
	y = ax + b

	parameters
	-------------
	X  : independent variable data (array form)
	Y  : dependent variable data (array form)

	return
	-------------
	a  : the slope of regression line
	b  : intercept
	r2  : coefficient of determination

	"""

	n = len(X)
	sum_XY = np.dot(X,Y); sum_X = sum(X); sum_Y = sum(Y); sum_X2 = sum(X**2); sum_Y2 = sum(Y**2)

	a = (n*sum_XY - sum_X*sum_Y) / (n*sum_X2 - (sum_X)**2)
	b = (sum_Y - a*sum_X) / n

	r2 = (n*sum_XY - sum_X*sum_Y)**2 / ((n*sum_X2 - sum_X**2) * (n*sum_Y2 - sum_Y**2))

	return a, b, r2