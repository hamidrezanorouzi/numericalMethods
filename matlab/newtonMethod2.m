% ---------------- Matlab function ---------------------------
% Numerical methods course, AUT
% website: www.cemf.ir
% newton's method to one find root of fx 
% DEFs:
%inputs:
%   fx: function
%   x0: initial guess of root
%   tol: relative tolerance 
%output:
%   x1: root
%   fx1: function value at root
%   iter: number of iterations
%   ea: approximate error 

function [x1, fx1, iter, ea] = newtonMethod2(fx, x0, tol)

if( nargin <3 )
    tol = 1.0e-5; %default relative tolerance 
end

%interval for numerical drivatives
dx = max(0.00001, tol);

nMax = 100; %maximum number of iterations
fx0 = fx(x0);

for iter = 1:nMax
    
    fx0dx = fx(x0+dx);
    dfx0 = (fx0dx - fx0)/dx; %forward derivatives
    
    x1 = x0 - fx0/dfx0;
    fx1 = fx(x1);
    
    if( abs(fx1) < tol)
        ea = abs(fx1);
        return;
    end
    
    if( abs(x1)>1.0e-15)
        ea = abs((x1-x0)/x1);
    else
        ea = abs(x1);
    end
    
    if( ea < tol)
       return;
    end
    x0 = x1;
    fx0 = fx1;
end

% no solution 
error('newtonMethod2 did not converged in %d iteration', nMax)

end

