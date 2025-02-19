# Matrix Functions

A Python package for **numerically computing matrix functions**, with a particular focus on **difference equations** and **analytic continuation of matrices**.

## Features

- **General Matrix Functions** – Computes **any function** of a matrix.
- **Numerical Computation** – Focuses on **floating-point approximations** rather than symbolic computation.
- **Difference Equations** – Provides tools for solving **recurrence relations** using matrix function techniques.
- **Analytic Continuation** – Enables non-integer shifts in difference equations using analytic continuation.
- **Mathematical Documentation** – Each release includes a **PDF document** explaining the mathematical foundations.
- **Pure Python** – No compilation required, making installation simple and cross-platform.

## Installation

The package is available on PyPI and can be installed with:

```bash
pip install matrixfuncs
```

###

## Usage

The Cayley-Hamilton theorem states that every square matrix is a root of its own characteristic polynomial. From this it follows that $A^n$ is a linear combination of $A^0,\ A,\ A^2,\ \dots,\ A^{n-1}$ and therefore every polynomial in A is such a linear combination:

$A^m = \sum_{k=0}^{n-1} \alpha_{mk} A^k$

It turns out that $\alpha_{mk}$ only depends on the eigenvalues $\lambda_1,\ \dots,\ \lambda_n$. Hence, every matrix function can be expressed in in such a way if the function is analytic in the eigenvalues:

$f(A) = \varphi_{ij}^{(k)} A^i f^{(k)}(\lambda_j)$


### Computing the Sine Function Using Matrix Recurrence Relations

The sine function satisfies the recurrence relation:

$sin(x + a) = 2\cos(a) sin(x) - sin(x - a)$

This allows us to express the sine of a shifted angle using matrix multiplication:

$
\begin{bmatrix} sin(x + a) \\ sin(x) \end{bmatrix} =
\begin{bmatrix} 2\cos(a) & -1 \\ 1 & 0 \end{bmatrix}
\begin{bmatrix} sin(x) \\ sin(x - a) \end{bmatrix}
$



Using **matrix exponentiation**, we can compute \(sin(x + na)\) efficiently.

#### Python Implementation

```python
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
fig = plt.figure(figsize=(8, 5))
plt.plot(ts, f_M, 'b-', label='Continuation of sampled function')
plt.plot(unit*np.arange(2*np.pi/unit), np.sin(unit*np.arange(2*np.pi/unit)), 'ro', label=f'Sampled at step size {unit}')
plt.xlabel('Time ($t$)')
plt.ylabel('$\\sin(t)$')
plt.title(f'Smooth Continuaton of the Sine function with Step Size {unit}')
plt.legend()
plt.grid(True)
plt.show(fig)
fig.savefig('sinFromMFunc.png', dpi=200)
```
![Output of plt.plot(ts, sinFromMFunc)](https://raw.githubusercontent.com/nextdorf/matrix-functions/main/sinFromMFunc.png?raw=true)


## Documentation

- A **PDF document** explaining the mathematical background is included with each release.
- The package includes **in-code documentation** but no separate programming guide at this time.

## Benchmarks & Performance

Currently, the package is focused on **numerical accuracy** rather than high-performance computing. Future updates will prioritize:

- **Stability improvements** for larger matrices.
- **Optimized numerical methods** to reduce floating-point errors.
- **Support for large-scale computations**.

## Future Plans

This package is **not yet feature-complete**. Future improvements will focus on:

- Enhancing **numerical stability**.
- Expanding support for **larger matrices**.
- General optimizations and performance improvements.

## Contributing

Contributions are welcome! If you'd like to help improve this package:

1. Fork the repository.
2. Create a new branch for your changes.
3. Submit a pull request with a clear description of your updates.

Areas where contributions could be particularly helpful:

- **Performance optimization**.
- **Expanding function support**.
- **Adding better documentation**.

For feature suggestions or bug reports, please open an issue on GitHub.

## License

This project is licensed under the **LGPL-3 License**. See the [LICENSE](LICENSE) file for details.