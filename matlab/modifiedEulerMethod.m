% ---------------- Matlab function ---------------------------
% Numerical methods course, AUT
% website: www.cemf.ir
% Modified Euler method to solve an ivp ODE
% DEFs:
%inputs:
%   dydt: function to be integrated
%   tspan: [ti, tf] where ti and tf = initial and final values of independent variable
%   y0: initial value of dependent variable
%   n: number of steps
%   p1,p2,... = additional parameters used by dydt
%output:
%   t: vector of independent variable
%   y: vector of solution for dependent variable

function [t,y] = modifiedEulerMethod(dydt, tspan, y0, n, varargin)

    if (nargin<4)
        error('at least 4 input arguments required');
    end
    ti = tspan(1);tf = tspan(2);  
    
    if (~(tf>ti))
        error('upper limit must be greater than lower');
    end
    h = (tf-ti)/n;
    y = zeros(n+1,1);
    t = zeros(n+1,1);
    tol = 1.0e-4;
    %initial condition
    t(1) = ti;
    y(1) = y0;
    
    for i = 1:n %implement Heun's method
        t(i+1) = t(i) + h; fxi = dydt( t(i), y(i), varargin{:} );
        
        %predictor
        yold = y(i)+h*fxi;
        
        converged = false;
        for iter = 1:1000
           fxi1 =  dydt( t(i+1), yold, varargin{:} );
           ynew = y(i)+ h*(fxi+fxi1)/2.0;
           if( abs((ynew-yold)/max(abs(ynew),tol)) < tol) %avoids divided by zero
               converged = true;
               y(i+1) = ynew;
               break;
           else
               yold = ynew;
           end
        end
        
        if(~converged)
            error('solution of predictor-corrector did not converged in 1000 iterations');
        end
    end
    
end