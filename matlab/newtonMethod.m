% ---------------- Matlab function ---------------------------
% Numerical methods course, AUT
% website: www.cemf.ir
% newton's method to one find root of fx 
% DEFs:
%inputs:
%   fx: function
%   dfx: derivative of function 
%   x0: initial guess of root
%   tol: relative tolerance 
%output:
%   x1: root
%   fx1: function value at root
%   iter: number of iterations
%   ea: approximate error 

function [x1, fx1, iter, ea] = newtonMethod(fx, dfx, x0, tol)

if( nargin <4 )
    tol = 1.0e-5; %default relative tolerance 
end


nMax = 100; %maximum number of iterations


for iter = 1:nMax
    x1 = x0 - fx(x0)/dfx(x0);
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
end

% no solution 
error('newtonMethod did not converged in %d iteration', nMax)

end

