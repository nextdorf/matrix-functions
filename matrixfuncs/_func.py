from typing import Callable
from ._krylov import KrylovSpace
import numpy as np


class MFunc:
  '''A helper class to calculate a matrix function using the krylov space

  Attributes
  ----------
  coeffs : array-like
    The coefficients used to calculate the matrix function.
    If `appliedFunc` is False then `coeffs_(...k)f(ev_k) = T_(...ij)f(A)_ij`, otherwise `coeffs_(...) = T_(...ij)f(A)_ij`. For some tensor `T` and `ev` being the eigenvalues.
  eigvals : array-like
    The eigenvalues of the matrix
  appliedFunc : bool = False
    Specifies whether a function was already applied by contraction to the coefficients.

  Example
  -------
  \\# One can show that sin(x) = 2cos(a)sin(x-a) - sin(x-2a). So for `x=0` and `a=.2` we can define a matrix `M`
  such that [sin(x+na), sin(x+(n-1)a)] = matrix_power(M, n) @ [sin(x), sin(x-a)]:
  >>> import numpy as np
  >>> a = .2
  >>> M = np.array([[2*np.cos(a), -1], [1,0]])
  >>> for n in [0, 1, 10, 100, 1000, 10000]:
  >>>   sinNa = [1,0] @ np.linalg.matrix_power(M, n) @ np.sin([0,-a])
  >>>   error = (sinNa - np.sin(n*a))**2
  >>>   print('n:', error)
  >>> ...
  n: 0.0
  n: 0.0
  n: 4.930380657631324e-30
  n: 1.199857436841159e-27
  n: 9.247078067402441e-26
  n: 3.397767010759534e-24

  \\# In order to generalize the discrete `n` to a continous variable we have to replace `matrix_power(M, n)` with
  something like `exp(n*log(M))`. For calculating the function of a matrix we can use `MFunc`.
  >>> krylovM = KrylovSpace(M)
  >>> mfn = MFunc.fromKrylov(krylovM)
  >>> mfn0 = [1,0] @ mfn @ np.sin([0,-a])

  \\# Compare the numeric evaluation of sin and numpy's implementation in 1000 points between 0 and 2pi:
  >>> ts = np.linspace(0, 2*np.pi, 1000)
  >>> sinFromMFunc = mfn0(lambda evs: np.exp(np.outer(np.log(evs), ts/a))).coeffs
  >>> sints = np.sin(ts)
  >>> error = np.sum((sinFromMFunc - sints)**2)
  >>> print(error)
  1.1051029979680703e-26
  '''
  def __init__(self, coeffs: np.ndarray, eigvals: np.ndarray, appliedFunc: bool = False):
    self.eigvals = eigvals
    self.coeffs = coeffs
    self.appliedFunc = appliedFunc

  def __call__(self, fn: Callable, *args, **kwds):
    '''Applies the function fn to the coefficients and returns the resulting `MFunc`
    Apply like `coeffs_(...k)f(ev_k) = T_(...ij)f(A)_ij`

    fn : array-like -> array-like
      Should evaluate the eigenvalues such that ev_k is mapped to f(ev_k)_(k...).

    args, kwds
      Additional arguments are passed to fn.
    '''
    if self.appliedFunc:
      return None
    fnEigvals = fn(self.eigvals, *args, **kwds)
    appliedCoeffs = np.tensordot(self.coeffs, fnEigvals, ([-1], [0]))
    return MFunc(np.real_if_close(appliedCoeffs), self.eigvals, True)

  @staticmethod
  def fromKrylov(space: KrylovSpace):
    'Constructs `MFunc` from a `KrylovSpace` object'
    return MFunc(space.funcCoeffs.copy(), space.eigvals)

  def __matmul__coeffs(self, other):
    if self.appliedFunc:
      coeffs = np.tensordot(self.coeffs, other, ([-1], [0]))
      if isinstance(coeffs, type(NotImplemented)):
        return coeffs
    else:
      coeffs = np.tensordot(self.coeffs, other, ([-2], [0]))
      if isinstance(coeffs, type(NotImplemented)):
        return coeffs
      oldRank = len(self.coeffs.shape)
      newRank = len(coeffs.shape)
      if oldRank-1 != newRank:
        coeffs = np.moveaxis(coeffs, oldRank-2, -1)
    return coeffs
  def __matmul__(self, other):
    coeffs = self.__matmul__coeffs(other)
    return MFunc(coeffs, self.eigvals, self.appliedFunc) if not isinstance(coeffs, type(NotImplemented)) else coeffs
  def __rmatmul__(self, other):
    coeffs = np.tensordot(other, self.coeffs, ([-1], [0]))
    if isinstance(coeffs, type(NotImplemented)):
      return coeffs
    return MFunc(coeffs, self.eigvals, self.appliedFunc)
  def __imatmul__(self, other):
    coeffs = self.__matmul__coeffs(other)
    if isinstance(coeffs, type(NotImplemented)):
      return coeffs
    self.coeffs = coeffs
    return self

  def __repr__(self) -> str:
    return 'MFunc'+repr(dict(Eigvals = self.eigvals, Coeffs=self.coeffs))


