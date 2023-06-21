function x = gauss_seidel(A, B, X0, tol, max_iter)

% A - coefficient matrix
% X0 - initial guess for solution vector
% tol - convergence tolerance
% max_iter - maximum number of iterations
% Output: x - approximate solution vector
% nr is the number of the rows of A.
% nc is the number of the columns of B.

% Check if A is square and B is a column vector
[nr, nc] = size(A);
if nr ~= nc
    error('A must be a square matrix')
end
if size(B, 2) ~= 1
    error('B must be a column vector')
end

% Initialize x and error
x = X0;
err = inf(nr, 1);

% Perform Gauss-Seidel iteration
iter = 0;
while any(abs(err) > tol) && iter < max_iter
    iter = iter + 1;
    for i = 1:nr
        % Compute the sum of A(i,j)*x(j) for j ~= i
        sum = 0;
        for j = [1:i-1 i+1:nr]
            sum = sum + A(i,j)*x(j);
        end
        % Update x(i) using the Gauss-Seidel formula
        x(i) = (B(i) - sum)/A(i,i);
    end
    % Compute the error
    err = B - A*x;
end

% Display the number of iterations and the error norm
fprintf('Number of iterations: %d\n', iter)
fprintf('Error norm: %g\n', norm(err))
end
