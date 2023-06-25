import numpy as np
import matplotlib.pyplot as plt

def r_squared(x, y):
    # Calculate the mean of y values
    y_mean = np.mean(x)

    # Calculate the total sum of squares
    total = np.sum((x - y_mean) ** 2)

    # Calculate the residual sum of squares
    residual = np.sum((x - y) ** 2)

    # Calculate the R-squared value
    rsquared = 1 - (residual / total)

    return rsquared

# Input data
x = np.array([1, 2, 3, 4, 5])
y = np.array([100, 85, 80, 70, 70])

# Fit a polynomial model of degree 3
mymodel = np.poly1d(np.polyfit(x, y, 3))

# Generate a line of values within the range of x
myline = np.linspace(x.min(), x.max(), 100)

# Calculate the R-squared value of the model fit
fit = r_squared(y, mymodel(x))

# Print the R-squared value
print("The R-squared is %f. The closer it is to one, the more accurate the model is." % fit)

# Plot the scatter plot of the data points
plt.scatter(x, y)

# Plot the polynomial model
plt.plot(myline, mymodel(myline))

# Display the plot
plt.show()

