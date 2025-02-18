import matrixfuncs as mf
import numpy as np
import matplotlib.pyplot as plt

# Define a transformation matrix based on a unit parameter
unit = 0.2
M = np.array([[2 * np.cos(unit), -1], [1, 0]])

# Generate time steps for evaluation
ts = np.linspace(0, 2 * np.pi, 1000)

# Define the function to be applied to the matrix
f = lambda x: x ** (ts / unit)

# Convert the matrix into a functional array representation
arr = mf.FArray.from_matrix(M)

# Define input vectors for left-hand side and right-hand side multiplications
v0_lhs = np.array([1, 0])  # Left-hand side vector
v0_rhs = np.sin([0, -unit])  # Right-hand side vector

# Compute the function applied to the matrix and evaluate it
vals = v0_lhs @ arr @ v0_rhs  # Compute matrix function application
f_M = vals(f)  # Evaluate the function over the time steps


# Plot the computed function values
plt.figure(figsize=(8, 5))
plt.plot(ts, f_M, 'b-', label='Continuation of sampled function')
plt.plot(unit*np.arange(2*np.pi/unit), np.sin(unit*np.arange(2*np.pi/unit)), 'ro', label=f'Sampled at step size {unit}')
plt.xlabel('Time ($t$)')
plt.ylabel('$\\sin(t)$')
plt.title(f'Smooth Continuaton of the Sine function with Step Size {unit}')
plt.legend()
plt.grid(True)
plt.show()
plt.savefig('sinFromMFunc.png', dpi=200)
