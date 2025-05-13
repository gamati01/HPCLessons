import sys
import numpy as np
import matplotlib.pyplot as plt

N = sys.argv[1]

# --- Parameters: adjust these based on your Fortran simulation ---
Nx = int(N)+2  # Number of grid points in x (should match your Fortran code)
Ny = int(N)+2  # Number of grid points in y

# --- Load the data ---
# The file 'laplace_solution.dat' is assumed to have three columns: x, y, u
data = np.loadtxt('laplace_solution.dat')

# Extract the columns
x = data[:, 0]
y = data[:, 1]
u = data[:, 2]

# Reshape the 1D arrays into 2D arrays.
# This assumes the data is written row by row (i.e. each row corresponds to a fixed y value).
X = x.reshape((Ny, Nx))
Y = y.reshape((Ny, Nx))
U = u.reshape((Ny, Nx))

EXACT = np.sin(np.pi*X)*np.exp(-np.pi*Y)

# --- Plot the solution ---
plt.figure(figsize=(8, 6))
contour = plt.contourf(X, Y, U, 12, cmap='viridis')
plt.contour(X, Y, EXACT, 12, colors='white')
plt.colorbar(contour, label='u')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Solution of the Laplace Equation')
plt.tight_layout()
plt.show()

