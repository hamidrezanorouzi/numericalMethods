% ---------------- Matlab function ---------------------------
% Numerical methods course, AUT
% website: www.cemf.ir
% newton's method to one find root of fx 
% DEFs:
%inputs:
%   fx: function
%   x0: a vector of initial guess
%   tol: relative tolerance 
%   alpha: under-relaxation
%output:
%   x1: solution of the set
%   fx1: functions values at root
%   iter: number of iterations
%   ea: approximate error 

function [x1, fx1, iter, ea] = newtonMethodSet2(fx, x0, tol, alpha)

    if( nargin <3 )
        tol = 1.0e-5; %default relative tolerance 
        alpha = 1;
    elseif( nargin <4)
        alpha = 1;
    end

    %interval for numerical drivatives
    dx = max(0.00001, tol);

    nMax = 100; %maximum number of iterations

    % start newton loop
    for iter = 1:nMax
               
        [jac , f] = jacobianMat(fx, x0, dx);
        
        % check if x0 is the solution of the set
        ea = dot(f,f);
        if ( ea < tol )
            fx1 = f;
            x1 = x0;
            return;
        end
        
        delta = jac\(-f);
        x1 = x0 + alpha*delta;
        x1
        ea = dot( delta, delta )/ max( max(x1), 1.0 ) ;
        if( ea  < tol)
            fx1 = fx(x1);           
            return;
        end
        
        x0 = x1;
    end
    
    % no solution 
    error('newtonMethodSet2 did not converged in %d iteration', nMax)

end


function [jac, fx0] = jacobianMat(fx, x, dx)

    n = length(x); %number of equations
    jac = zeros(n,n);
    fx0 = fx(x);
    
    for i=1:n
        holdXi = x(i);
        x(i) = x(i)+ dx;
        fx1 = fx(x);
        x(i) = holdXi;
        jac(:,i) = (fx1-fx0)./dx; %derivative with respect to xi
    end
    
    return;
end
