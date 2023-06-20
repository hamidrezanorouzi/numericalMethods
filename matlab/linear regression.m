% Sample data
x = [1, 2, 3, 4, 5];
y = [2, 4, 5, 4, 5];

% Calculate the mean of x and y
mean_x = mean(x);
mean_y = mean(y);

% Calculate the slope (b) and y-intercept (a)
numerator = 0;
denominator = 0;
for i = 1:length(x)
    numerator = numerator + (x(i) - mean_x) * (y(i) - mean_y);
    denominator = denominator + (x(i) - mean_x) ^ 2;
end
b = numerator / denominator;
a = mean_y - b * mean_x;

% Print the equation of the line
fprintf('y = %.2fx + %.2f\n', b, a);