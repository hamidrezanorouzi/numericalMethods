# Numerical methods course, AUT
# website: www.cemf.ir
import numpy as np

def gaussElimination(A, B, tol = 1.0e-12):
    """
    Gauss elimination method to solve a set of linear
    algebraic equations.
    inputs:
        A: nxn Coefficient matrix, should be from numpy
        B: nx1 vector of known values
    outputs:
        [1]: solution of the set
        [2]: determinant of matrix
    """
    #number of rows and columns
    nr,nc = A.shape;
    if nr!=nc:
        raise Exception('input matrix A is not square!')
    
    #number of equations in the set
    n,m = B.shape
    if n!= nr:
        raise Exception('A and B size mismatch')

    Aug = np.concatenate([A,B],1);
    detA = 1;
    
    #main loop
    for k in range(n-1):
        #first: partial pivoting
        pElement = abs(Aug[k,k])
        pRow = k;
        
        #locate maximum  element in the rows below the pElement
        for row in range(k+1,n):
            if abs(Aug[row,k]) > pElement:
                pElement = abs(Aug[row,k])
                pRow = row

        #interchanges the rows, if necessary
        if pRow != k:
            temp = np.array(Aug[k,:])
            Aug[k,:] = Aug[pRow,:]
            Aug[pRow,:] = temp
            detA = -detA  #change of sign
        
        if abs(Aug[k,k]) < tol:
            raise Exception('Singular matrix!')
        

        # making elements below the pivot element zero
        for m in range(k+1,n):
            Aug[m,:] -= Aug[m,k]/Aug[k,k] * Aug[k,:]
    
    #the last equation never checked
    if abs(Aug[n-1,n-1]) < tol:
        raise Exception('Singular matrix!')

    # back substitution
    X = np.zeros((n,1))
    X[n-1] = Aug[n-1,n]/Aug[n-1,n-1];
    detA = detA*Aug[n-1,n-1];
    
    for k in range(n-2, -1, -1):
        X[k] = (Aug[k,n] - np.dot(Aug[k,k+1:n],X[k+1:n]) )/Aug[k,k]
        detA = detA * Aug[k,k]

    return (X, detA);
