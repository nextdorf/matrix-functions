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
