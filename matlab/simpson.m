% Define the function to be integrated
f = @(x)(0.2 + 25*x - 200*x^2 + 675*x^3 - 900*x^4 + 400*x^5);

% Define the integration limits
% enter your lower limit here
a =1.4;

% enter your upper limit here
b =2.2;

% enter the number of subintervals (must be even)
n = 9;

% Call the simpson function to get the approximation
I = simpson(f,a,b,n);

% Display the result
fprintf('The approximate value of the integral is %f\n',I);

% Define the simpson function
function I = simpson(f,a,b,n)
    % Check if n is 3k
    if mod(n,3) ~= 0
        error('n must be even');
    end
    
    % Calculate the step size
    h = (b-a)/n;
    
    % Initialize the sum
    sum = 0;
    
    % Loop over the subintervals
    for i = 1:n/3
        % Get the endpoints and the midpoint of the subinterval
        x0 = a + 3*(i-1)*h;
        x1 = x0 + h;
        x2 = x1 + h;
        x3 = x2 + h;
        
        % Apply Simpson's 3/8 rule to the subinterval
        sum = sum + 3*h/8 * (f(x0) + 3*f(x1) + 3*f(x2) + f(x3));
    end
    
    % Return the final approximation
    I = sum;
end