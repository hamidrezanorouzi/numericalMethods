% ---------------- Matlab function ---------------------------
% Numerical methods course, AUT
% website: www.cemf.ir
% Root finding based on the bisection method in interval (a b)
% DEFs:
% bisection(fx, a, b): finds the root of fx in (a b).
% bisection(fx, a, b, tol): finds the root of fx in (a b) with
%   relative tolerance tol.
function [x0, fx0, iter, ea] = bisection(fx, a, b, tol)

if( nargin <4 )
    tol = 1.0e-5;
end

nMax = 1000; 

fa = fx(a);
fb = fx(b);
c_old = a;

%main loop
for iter=1:nMax
    
    c = (a+b)/2;
    fc = fx(c);
    if( abs(fc) < tol)
        x0 = c;
        fx0 = fc; ea = fc;
        return;
    end
        
    if( fc*fb < 0 )
        a = c;
        fa = fc;
    else
        b = c;
        fb = fc;
    end
    
    if( abs(c) > 1.0e-15 )
        ea = abs((c-c_old)/c);
    else
        ea = abs(c);
    end
    
    if( ea < tol)
        x0 = c;
        fx0 = fc;
        return;
    end
    c_old = c;
end

% no solution 
error('bisection did not converged in %d iteration', nMax)

end


