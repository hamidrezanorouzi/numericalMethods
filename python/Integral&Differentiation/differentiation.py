import numpy as np 


def CenteralWeights(n = 1 , h_order = 1):

	"""
	return weight of ceteral formula

	if n=1 and h_order = 2:
	f'(xi) = w[0]* f(xi+2) + w[1]*f(xi+1) + w[2]*f(xi-1) + w[3]*f(xi-2)
	these multiply must done by another function

	parameter
	--------------
	n : is order of differentiation is 1 to 4
	h_order : 2 / 4 ==> O(h2) / O(h4)

	Notes
    -------------
    Decreasing the step size too small can result in round-off error.
    instead use RombergDiff

	"""
	if n > 4: raise Exception("more than forth derivative is not implement")  
	if h_order not in [2,4] : raise Exception("only O(h**4) and O(h**2) is implement")

	numinator_weight = [[[-1, 0, 1], [1, -8, 0, 8, -1]],
				  [[1, -2, 1], [-1, 16, -30, 16, -1]],
				  [[-1, 2, 0, -2, 1], [-1, 8, -13, 0, 13, -8, 1]],
				  [[1, -4, 6, -4, 1],[-1, 12, -39, 56, -39, 12, -1]]
				  ]
	dominator_weight = [[[2], [12]], 
						[[1], [12]],
						[[2], [8]],
						[[1], [6]]
					   ]

	i = n-1; j = int(h_order/2) -1

	return np.array(numinator_weight[i][j]), np.array(dominator_weight[i][j])



def CenteralDifferentiation(f, x, h = 0.1, n = 1, h_order = 4):

	"""	
	estimate of n-th derivative of f(x) with centeral formula

	parameters
	-------------
	f  : function we want it's derivative
	x  : where we want derivative to calculate
	h  : step size
	n  : order of derivative (first to 4th derivative)
	h_order : 2 / 4 ==> O(h2) / O(h4) 

	return
	-------------
	derivative estimation

	"""

	num, dom = CenteralWeights(n, h_order)

	N = len(num) # number of function must evaluate
	a = x - (N-1)/2 * h;  b = x + (N-1)/2 * h
	x = np.linspace(a, b, N)

	diff = sum(num * f(x) / (dom * h**n))

	return diff
