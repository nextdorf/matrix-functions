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

$$
A^m = \sum_{k=0}^{n-1} \alpha_{mk} A^k
$$

It turns out that $\alpha_{mk}$ only depends on the eigenvalues $\lambda_1,\ \dots,\ \lambda_n$. Hence, every matrix function can be expressed in in such a way if the function is analytic in the eigenvalues:

$$
f(A) = \varphi_{ij}^{(k)} A^i f^{(k)}(\lambda_j)
$$


### Computing the Sine Function Using Matrix Recurrence Relations

The sine function satisfies the recurrence relation:

$$
sin(x + a) = 2\cos(a) sin(x) - sin(x - a)
$$

This allows us to express the sine of a shifted angle using matrix multiplication:

$$
\begin{bmatrix} sin(x + a) \\ sin(x) \end{bmatrix} =
\begin{bmatrix} 2\cos(a) & -1 \\ 1 & 0 \end{bmatrix}
\begin{bmatrix} sin(x) \\ sin(x - a) \end{bmatrix}
$$



Using **matrix exponentiation**, we can compute \(sin(x + na)\) efficiently.

#### Python Implementation

```python
import numpy as np
from matrixfuncs import *
import matplotlib.pyplot as plt

# Define the parameter 'a'
a = 0.2

# Construct the matrix M
M = np.array([
  [2 * np.cos(a), -1],
  [1, 0]
])

# Compute and compare sin(na) for n=0, 1, 10, 100, 1000, 10000
for n in [0, 1, 10, 100, 1000, 10000]:
  # sin(n a) using numpy
  sinNa = [1,0] @ np.linalg.matrix_power(M, n) @ np.sin([0,-a])

  # squared error between directly computing np.sin and using the recurrence relation via M^n
  error = (sinNa - np.sin(n*a))**2

  print('n:', error)

# n: 0.0
# n: 0.0
# n: 4.930380657631324e-30
# n: 1.199857436841159e-27
# n: 9.247078067402441e-26
# n: 3.397767010759534e-24

# In order to generalize the discrete n to a continous variable we have to replace
# matrix_power(M, n) with something like exp(n*log(M)). For calculating the function of a matrix
# we can use MFunc.
krylovM = KrylovSpace(M)
mfn = MFunc.fromKrylov(krylovM)
mfn0 = [1,0] @ mfn @ np.sin([0,-a])

# Compare the numeric evaluation of sin to numpy's implementation in 1000 points between 0 and 2pi:
ts = np.linspace(0, 2*np.pi, 1000)
sinFromMFunc = mfn0(lambda evs: np.exp(np.outer(np.log(evs), ts/a))).coeffs
sints = np.sin(ts)
error = np.sum((sinFromMFunc - sints)**2)
print(error)
# 1.1051029979680703e-26

plt.plot(ts, sinFromMFunc)
plt.savefig('sinFromMFunc.png', dpi=200)
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

This project is licensed under the **GPL-3 License**. See the [LICENSE](LICENSE) file for details.