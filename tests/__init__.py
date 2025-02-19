from matrixfuncs import *
from matrixfuncs.utils import *
import numpy as np
import pytest

def direct_prod(xs0, *xss):
  if not hasattr(xs0, '__iter__'):
    xs0 = (xs0, )
  if len(xss) == 0:
    for x in xs0:
      yield (x,)
  else:
    inner = list(direct_prod(*xss))
    for x in xs0:
      for y in inner:
        yield (x,) + y

def matrix_from_mult(mult: Multiplicity):
  t=np.arange(mult.dim)
  M=np.zeros(t.shape*2)
  M[t,t] = [ev for _,ev,_,_ in mult.ev_iter('full')]
  off_diag0 = [1.*(np.arange(alg) < alg-geom) for alg, geom in zip(mult.algebraic, mult.geometric)]
  off_diag = np.concat(off_diag0)[:-1]
  M[t[:-1],t[1:]] = off_diag
  return M

def rnd_matrix_from_mult(mult: Multiplicity):
  M0 = matrix_from_mult(mult)
  Q, Q_inv = None, None
  while Q_inv is None:
    try:
      if np.iscomplexobj(M0):
        Q = np.random.randn(*M0.shape, 2) @ [1, 1j]
      else:
        Q = np.random.randn(*M0.shape)
      Q /= np.linalg.det(Q)
      Q_inv = np.linalg.inv(Q)
    except:
      continue
    inv_err = err(Q_inv@Q, np.eye(M0.shape[0]))
    if inv_err > 1e-9:
      pytest.warns(f'Too big error in rnd_matrix_from_mult ({inv_err}), try again...')
      Q_inv = None
  M = Q @ M0 @ Q_inv
  return M


