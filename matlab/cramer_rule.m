function x=cramer_rule(A,B)

% this function is used to solve linear system of equation.
% the general form is Ax=B.
% x is the matrix of unknown variables.
% A is the coefficient matrix of x.
% nr is the number of the rows of A.
% nc is the number of the columns of B.

[nr,nc]=size(A);  % calculating the size of A.
n=length(B);      % calculating the length of B.
D=det(A);         % calculating the determinant of A.

% conditions 
% if these conditions are uncertaint,cramer's rule can't solve the system.

if nr~=nc         % checking if the matrix is square.
    fprintf('error!coefficient matrix is not square')
end

if n~=nr          % condition of the possibblity of multiplication.
    fprintf('error!can not multiple coefficient matrix and vector')
end

if D==0           % checking the determinant of A if it is not zero.
    fprintf('error!can not be answered by cramer rule')
end

for i=1:n
    Z=A;          % an assistant matrix.
    Z(:,i)=B;     % changing B with the columns of A.
    d=det(Z);
    x(i)=d/D;
end
end