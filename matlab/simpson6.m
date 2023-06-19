% Simpson's 1/3 rule function
function y = simpson6(f,a,b,n)
% f is the function to be integrated
% a and b are the lower and upper limits of integration
% n is the number of sub-intervals (must be even)
% y is the approximate value of the integral

h = (b-a)/n; % calculate the step size
x = a:h:b; % create an array of x values
y = f(a) + f(b); % initialize the sum with the first and last terms
for i = 2:n % loop through the odd indices from 2 to n
    if mod(i,2) == 0 % check if i is even
        y = y + 4*f(x(i)); % add four times the function value at x(i)
    else % otherwise i is odd
        y = y + 2*f(x(i)); % add two times the function value at x(i)
    end
end
y = y * h/3; % multiply the sum by h/3 to get the final answer
end