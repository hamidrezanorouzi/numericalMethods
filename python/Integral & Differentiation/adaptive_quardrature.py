from simpson import *


def Quadrature(f, span, n = 2, err = 1.49e-8):

	"""
	both romber and simpson's rule use equally spaced
	points. This perspective ignores the fact that 
	many functions have regions of high variability 
	along with other section where change is gradual

	in Adapive Quadrature we calculate integral with
	simpson 1_3 with n and n+1 segment and if the 
	difference between those calculation is large 
	we need to add segment 

	parameters 
	-------------

	f : function in integral
	span : [a, b] the span of integral
	n : initial number of segment
	err : the criterion for return |I2 - I1| <= err

	return
	-------------

	estimate of integral

	"""
    I1 = simpson1_3(f, span, n)
    I2 = simpson1_3(f, span, n+2)

    if abs(I2 - I1) <= err:

        return I2 + 1/15 * (I2-I1)

    else:

        return quadrature(f, span , n+2)




