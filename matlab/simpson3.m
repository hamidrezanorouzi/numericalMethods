% Define the function to be integrated
f = @(x)(0.2 + 25*x - 200*x^2 + 675*x^3 - 900*x^4+400*x^5);

% Define the integration limits
a = 3;
b = 6;

% Define the number of subintervals (must be even)
n = 6;

% Call the simpson function to get the approximation
I = simpson(f,a,b,n);

% Display the result
fprintf('The approximate value of the integral is %f\n',I);

% Define the simpson function
function I = simpson(f,a,b,n)
    % Check if n is even
    if mod(n,2) ~= 0
        error('n must be even');
    end
    
    % Calculate the step size
    h = (b-a)/n;
    
    % Initialize the sum
    sum = 0;
    
    % Loop over the subintervals
    for i = 1:n/2
        % Get the endpoints and the midpoint of the subinterval
        x0 = a + 2*(i-1)*h;
        x1 = x0 + h;
        x2 = x1 + h;
        
        % Apply Simpson's 1/3 rule to the subinterval
        sum = sum + h/3 * (f(x0) + 4*f(x1) + f(x2));
    end
    
    % Return the final approximation
    I = sum;
end

