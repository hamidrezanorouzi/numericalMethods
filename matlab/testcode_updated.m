% Definition of coefficient matrix A and constant vector b
A=[2 3 -1 4;4 -1 0 6;1 -3 2 -1;3 4 -2 1]
b=[10;4;-5;6]
% Gauss_Jordan_updated
x = Gauss_Jordan_updated(A, b);

x_expected = A\b;

if isequal(x, x_expected)
    disp('code is correct');
else
     disp(' code is  correct ');
end
