import numpy as np
import scipy
from . import MSpace, apply_fn

class FArray:
  '''Represents a functional array associated with a matrix space.
  
  This class enables efficient function application to matrices using precomputed function coefficients.
  
  Attributes
  ----------
  space : MSpace
      The underlying matrix space associated with the functional array.
  mult_space : str
      Defines the multiplication space ('min' by default).
  bare : np.ndarray
      The function coefficients tensor. If not provided, it is copied from `space.f_coeffs`.
  '''
  __array_priority__ = 1
  def __init__(self, space: MSpace, mult_space='min', bare=None, stay_lifted=True):
    '''
    Initialize an FArray instance.
    
    Parameters
    ----------
    space : MSpace
        The matrix space associated with this functional array.
    mult_space : str, optional
        The multiplication space, default is 'min'.
    bare : np.ndarray, optional
        Precomputed function coefficients. Defaults to `space.f_coeffs` if not provided.
    stay_lifted : bool, optional
        Whether the object remains lifted when applying functions.
    '''
    self.space = space
    self.mult_space=mult_space
    self.bare = np.copy(space.f_coeffs) if bare is None else bare
    self.__stay_lifted = stay_lifted

  @staticmethod
  def from_matrix(M: np.ndarray, **kwargs):
    '''
    Create an FArray instance from a matrix.
    
    Parameters
    ----------
    M : np.ndarray
        The input matrix.
    **kwargs : dict
        Optional keyword arguments, including:
        - eigvals: Precomputed eigenvalues.
        - mult_space: Multiplication space setting.
        - bare: Precomputed function coefficients.
        - stay_lifted: Whether the instance should stay lifted.
    
    Returns
    -------
    FArray
        An instance of the functional array.
    '''
    _kwargs = {k: kwargs[k] for k in 'eigvals'.split() if k in kwargs}
    space = MSpace(M, **_kwargs)
    _kwargs = {k: kwargs[k] for k in 'mult_space bare stay_lifted'.split() if k in kwargs}
    return FArray(space, **_kwargs)

  def __call__(self, f, *dfs, gen_df=None, stay_lifted=None, real_if_close='auto', **kwargs):
    if isinstance(real_if_close, str):
      if real_if_close == 'auto':
        real_if_close = not np.iscomplexobj(self.space.M)
      else:
        raise ValueError(f'Unsupported value for real_if_close: {real_if_close}')

    ret = apply_fn(M=None, f=f, dfs=dfs, gen_df=gen_df, eigvals=self.space.multiplicity, coeffs=self.bare, mult_space=self.mult_space, real_if_close=real_if_close, **kwargs)

    ret_shape = np.shape(ret)
    is_squre_matrix = len(ret_shape)==2 and ret_shape[0] == ret_shape[1]
    if stay_lifted is None:
      stay_lifted = self.__stay_lifted and is_squre_matrix
    elif stay_lifted == True:
      stay_lifted = is_squre_matrix
    else:
      stay_lifted = False
    
    if stay_lifted:
      f_space = MSpace(ret)
      arr = FArray(f_space, mult_space=self.mult_space)
      return arr
    else:
      return ret

  @property
  def cs_shape(self):
    return np.shape(self.bare)[1:]
  @property
  def f_shape(self):
    return np.shape(self.bare)[0]
  @property
  def shape(self):
    return np.shape(self.bare)


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
    if not bool(self.cs_shape):
      return NotImplemented
    cs = np.tensordot(self.bare, other, 1)
    if cs is NotImplemented:
      return NotImplemented
    ret = FArray(self.space, mult_space=self.mult_space, bare=cs, stay_lifted=False)
    return ret
  def __rmatmul__(self, other):
    if not bool(self.cs_shape):
      return NotImplemented
    cs = np.tensordot(other, self.bare, ((-1,), (1,)))
    if cs is NotImplemented:
      return NotImplemented
    cs = np.moveaxis(cs, -len(self.cs_shape), 0)
    ret = FArray(self.space, mult_space=self.mult_space, bare=cs, stay_lifted=False)
    return ret
  def __imatmul__(self, other):
    if not bool(self.cs_shape):
      return NotImplemented
    cs = np.tensordot(self.bare, other, 1)
    if cs is NotImplemented:
      return NotImplemented
    self.bare = cs
    self.__stay_lifted = False
    return self

  def __repr__(self):
    ret = repr(self.bare)
    main_idx = ret.find('(')
    return f'FArray{ret[main_idx:]}'
