import numpy as np

def lu_decomposition(A):
    n = len(A)
    L = np.eye(n)
    U = np.zeros((n, n))

    for k in range(n):
        U[k, k:] = A[k, k:] - L[k, :k] @ U[:k, k:]
        L[k+1:, k] = (A[k+1:, k] - L[k+1:, :k] @ U[:k, k]) / U[k, k]

    return L, U

# Define the matrix
A = np.array([3, -0.1 ,-0.2], [0.1, 7, -0.3], [0.3, -0.2, 10]])

# Perform LU decomposition
L, U = lu_decomposition(A)

print("L (lower triangular matrix):\n", L)
print("U (upper triangular matrix):\n", U)
