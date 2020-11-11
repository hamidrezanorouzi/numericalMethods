% ---------------- Matlab function ---------------------------
% Numerical methods course, AUT
% website: www.cemf.ir
% Gauss elimination with partial pivoting to find the
% solution to AX=B 
% DEFs:
%inputs:
%   A: n x n coefficient matrix
%   B: n x 1 vector of known values
%output:
%   X: solution of the set
%   detA: determinant of A
function [X, detA] = gaussElimination(A,B, tol)

    %number of rows and columns    
    [nr, nc] = size(A); 
    
    if( nr~= nc)
        error('Coefficient matrix is not square!');
    end
    
    %number of equations in the set
    n = length(B);
    if( n~= nr)
        error( 'A and B size mismatch');
    end
    
    if( nargin<3)
        tol = 1.0e-12;
    end
    
    Aug = [A B]; %Augmented matrix
    detA = 1;
    
    %main loop for forward elimination 
    for k = 1:n-1
        %first: partial pivoting
        pElement = abs(Aug(k,k));
        pRow = k;
        %locate maximum element in the rows below the pElement
        for row = k+1:n
            if(abs(Aug(row,k)) > pElement)
                pElement = abs( Aug(row,k));
                pRow = row;
            end
        end
        
        %interchanges the rows, if necessary
        if( pRow ~= k)
            temp = Aug(k,:);
            Aug(k,:) = Aug(pRow,:);
            Aug(pRow,:) = temp;
            detA = - detA; %change of sign
        end
        
        if(abs(Aug(k,k)) < tol )
            error('Singular matrix!');
        end
        % making elements below the pivot element zero
        for m = k+1:n
            Aug(m,:) = Aug(m,:) - (Aug(m,k)/Aug(k,k)*Aug(k,:));
        end
        
    end
    
    %the last equation never checked
    if(abs(Aug(n,n)) < tol )
        error('Singular matrix!');
    end
    
    X = zeros(n,1);
    
    % back substitution
    X(n) = Aug(n,n+1)/Aug(n,n);
    detA = detA*Aug(n,n);
    
    for k=n-1:-1:1
        X(k) = (Aug(k,n+1) - dot(Aug(k,k+1:n),X(k+1:n)))/Aug(k,k);
        detA = detA * Aug(k,k);
    end
    
end 