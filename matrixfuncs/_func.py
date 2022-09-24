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
  '''
  def __init__(self, coeffs: np.ndarray, eigvals: np.ndarray, appliedFunc: bool = False):
    self.eigvals = eigvals
    self.coeffs = coeffs
    self.appliedFunc = appliedFunc

  def __call__(self, fn: Callable, *args, **kwds):
    '''Applies the function fn to the coefficients 
    '''
    if self.appliedFunc:
      return None
    fnEigvals = fn(self.eigvals, *args, **kwds)
    appliedCoeffs = np.tensordot(self.coeffs, fnEigvals, ([-1], [0]))
    return MFunc(np.real_if_close(appliedCoeffs), self.eigvals, True)

  @staticmethod
  def fromKrylov(space: KrylovSpace):
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


