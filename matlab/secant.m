% ---------------- Matlab function ---------------------------
% Numerical methods course, AUT
% website: www.cemf.ir
% Secant method to one find root of fx in [x0,x1]
% DEFs:
%inputs:
%   fx: function
%   x0: start point of search interval
%   x1: end point of search interval
%output:
%   x2: root
%   fx2: function value at root
%   iter: number of iterations
%   ea: approximate error 

function [x2, fx2, iter, ea] = secant( fx, x0, x1, tol)

    %inputs
    if( nargin <4 )
        tol = 1.0e-5; %default relative tolerance 
    end


    nMax = 1000; %maximum number of iterations
    x2_old = x0;

    fx0 = fx(x0);
    fx1 = fx(x1);

    for iter = 1:nMax
        x2 = x1 - (x1-x0)*fx1/(fx1-fx0);
        fx2 = fx(x2);

        if( abs(fx2) < tol) 
            ea = fx2;
            return;
        end
        if( fx2*fx1 < 0 )
            x0 = x2;
            fx0 = fx2;
        else
            x1 = x2;
            fx1 = fx2;
        end

        if( abs(x2) > 1.0e-15 )
            ea = abs((x2-x2_old)/x2);
        else
            ea = abs(x2);
        end

        if( ea < tol) 
            return;
        end

        x2_old = x2;

    end


    % no solution 
    error('secant method did not converged in %d iteration', nMax)

end