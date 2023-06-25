% Define the data points
x = [1, 2, 3, 4, 5];
y = [2, 4, 6, 8, 10];

% Define the degree of the polynomial fit
n = 2;

% Initialize the coefficient vector
p = zeros(n+1, 1);

% Calculate the polynomial fit using least squares method
for i = 0:n
    j = 0;
    s = 0;
    while j <= n
        if i+j == 0
            s = s + length(x);
        else
            s = s + sum(x.^i .* x.^j);
        end
        j = j + 1;
    end
    numerator = 0;
    for k = 1:length(x)
        if i == 0
            numerator = numerator + y(k);
        elseif i == 1
            numerator = numerator + y(k)*x(k);
        else
            numerator = numerator + y(k)*x(k)^i;
        end
    end
    p(i+1) = numerator / s;
end

% Plot the data and the polynomial fit
plot(x, y, 'o');
hold on;
plot(x, polyval(p, x));
xlabel('x');
ylabel('y');
legend('Data', 'Polynomial Fit');