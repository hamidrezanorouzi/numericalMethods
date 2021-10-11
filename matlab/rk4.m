% ---------------- Matlab function ---------------------------
% Numerical methods course, AUT
% website: www.cemf.ir
% Fourth order rk method to solve a set of ivp ODEs
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

function [t,y] = rk4(dydt, tspan, y0, n, varargin)

    if (nargin<4)
        error('at least 4 input arguments required');
    end
    ti = tspan(1);tf = tspan(2);  
    
    if (~(tf>ti))
        error('upper limit must be greater than lower');
    end
    nEq = length(y0);
    h = (tf-ti)/n;
    y = zeros(nEq,n+1);
    t = zeros(1, n+1);
    
    %initial condition
    t(1) = ti;
    y(:,1) = y0;
    
    for i = 1:n %implement rk method
        t(i+1) = t(i) + h;
        k1 = h * dydt( t(i), y(:,i), varargin{:});
        k2 = h * dydt( t(i)+0.5*h, y(:,i)+0.5*k1 , varargin{:});
        k3 = h * dydt( t(i)+0.5*h, y(:,i)+0.5*k2 , varargin{:});
        k4 = h * dydt( t(i)+ h, y(:,i)+ k3 , varargin{:});
        y(:,i+1) = y(:,i) + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
    end
    
end
