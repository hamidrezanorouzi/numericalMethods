% ---------------- Matlab function ---------------------------
% Numerical methods course, AUT
% website: www.cemf.ir
% Tomas method to solve a tri-diagonal set of linear equations 
% DEFs:
%inputs:
%   a: lower diagonal elements [1:n]
%   b: diagonal elements [1:n]
%   c: upper diagonal elements [1:n]
%   d: known coefficients [1:n]
%output:
%   X: solution of the set

function X = triDiagonal(a,b, c, d)
        
    %number of equations in the set
    n = length(b);
    
    cc = zeros(n,1);
    dd = zeros(n,1);
    
    %forward sweeping
    cc(1) = c(1)/b(1);
    for i=2:n-1
        cc(i) = c(i)/(b(i)-a(i)*cc(i-1));
    end
    
    dd(1) = d(1)/b(1);    
    for i=2:n
        dd(i) = (d(i)-a(i)*dd(i-1))/(b(i)-a(i)*cc(i-1));
    end
        
    X = zeros(n,1);
    % back substitution
    X(n) = dd(n);
    for i=n-1:-1:1
        X(i) = dd(i) - cc(i)*X(i+1);
    end
    
end 