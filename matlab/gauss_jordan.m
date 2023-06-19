function X = gauss_jordan(A,B)

% this function is used to solve linear system of equation.
% the general form is AX=B then converted to IX=C.
% X is the matrix of unknown variables.
% A is the coefficient matrix of X.
% nr is the number of the rows of A.
% nc is the number of the columns of B.

[nr,nc]=size(A);  % calculating the size of A.
n=length(B);      % calculating the length of B.

if diag(A)==0     % the diagonal array must be non zero.
    A=flipud(A);  % changing the rows of A.
    B=flipud(B);  % changing the rows of B.
end
        
D=det(A);         % calculating the determinant of A.

if nr~=nc         % checking if the matrix is square.
    fprintf('error!coefficient matrix is not square')
end
if D==0           % checking the determinant of A if is not zero.
    fprintf('error!can not be answered')
end
    
if n~=nr          % condition of the possibblity of multiplication.
    fprintf('error!can not multiple coefficient matrix and vector')
end

C=[A,B];          % sticking A to B by the columns

% start pivoting
    
for j=1:n           % for columns
    C(j,:)=C(j,:)/C(j,j);
    for i=1:n       % for rows
        if i~=j
            C(i,:)=C(i,:)-C(i,j)*C(j,:);
        end
    end
end

X=C(:,n+1);
end