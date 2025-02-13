import numpy as np
from collections import namedtuple
import numdifftools as nd
import scipy

class Multiplicity(namedtuple('Multiplicity', 'eigvals algebraic geometric'.split())):
  def __new__(cls, eigvals: np.ndarray, algebraic: np.ndarray, geometric: np.ndarray):
    return super(Multiplicity, cls).__new__(cls, eigvals, algebraic, geometric)

  def map(self, f, *dfs, gen_df=None, **kwargs):
    diff_count = np.max(self.algebraic)
    fs = [f] + list(dfs[:diff_count-1])
    if len(fs) < diff_count:
      if gen_df is None:
        k0 = len(fs)-1
        gen_df = lambda i, **_: nd.Derivative(fs[k0], n=i-k0, order=2*1)
      for i in range(len(fs), diff_count):
        fs.append(gen_df(i))
    ret = [fs[k](ev) for ev, alg in zip(self.eigvals, self.algebraic) for k in range(alg)]
    return np.array(ret)


def matrix_power_series(M: np.ndarray, stop: int):
  shape = np.shape(M)
  assert len(shape)==2 and shape[0]==shape[1], 'M is not a square matrix'
  dim = shape[0]
  ret = [np.eye(dim, dtype=M.dtype), M]
  for i in range(2, stop):
    k1 = i // 2
    k2 = i - k1
    ret.append(ret[k1] @ ret[k2])
  return np.array(ret)[:max(stop, 0)]

def eigval_multiplicity(M: np.ndarray, eigvals: np.ndarray | None =None, zero_thrsh = 1e-15, rel_eq_thrsh = 1e-8):
  if eigvals is None:
    eigvals = np.linalg.eigvals(M)
  else:
    eigvals = np.array(eigvals)
  nEigvals = eigvals.shape[0]
  non_zero_eigvals = eigvals[abs(eigvals) > zero_thrsh]
  unique_eigvals = [0*eigvals[0]] if non_zero_eigvals.shape != eigvals.shape else []
  for ev in non_zero_eigvals:
    differs = abs(np.array(unique_eigvals)/ev - 1) > rel_eq_thrsh
    if differs.all():
      unique_eigvals.append(ev)
  unique_eigvals = np.array(unique_eigvals)

  alg_mult = [0]*len(unique_eigvals)
  for ev in eigvals:
    alg_mult[abs(unique_eigvals - ev).argmin()] += 1
  alg_mult = np.array(alg_mult)
  
  geom_mult = []
  for i, ev in enumerate(unique_eigvals):
    if len(unique_eigvals) == 1:
      ev_thrsh = np.inf
    else:
      other_eigvals = unique_eigvals[np.arange(len(unique_eigvals)) != i]
      ev_thrsh = abs(other_eigvals - ev).min() / 2
    eigvals_ker, eigvecs_ker = np.linalg.eig(M - ev*np.eye(nEigvals))
    inds_ker = abs(eigvals_ker) < ev_thrsh
    vecs_ker = eigvecs_ker[:, inds_ker]
    # vecs_ker_min = normed_basis(vecs_ker.T).T
    vecs_ker_min = scipy.linalg.orth(vecs_ker)
    geom_mult.append(np.shape(vecs_ker_min)[-1])
  geom_mult = np.array(geom_mult)

  # Ret = namedtuple('Multiplicity', 'eigvals algebraic geometric'.split())
  # ret = Ret(unique_eigvals, alg_mult, geom_mult)
  ret = Multiplicity(unique_eigvals, alg_mult, geom_mult)
  return ret

def b_matrix(multiplicity):
  evs, alg, _ = multiplicity
  nEigvals = np.sum(alg)
  ret = np.zeros((nEigvals, nEigvals), dtype=evs.dtype)
  max_order = np.max(alg)
  inds0 = np.array([sum(alg[:i]) for i in range(len(alg))])
  val0 = np.concat([evs.reshape((-1, 1))**np.arange(nEigvals-1, 0, -1), np.ones((len(evs), 1), dtype=evs.dtype)], axis=1)
  for order in range(1, max_order+1):
    if order == 1:
      inds, val = inds0, val0
    else:
      evs_selection = alg >= order
      inds = (inds0 + order - 1)[evs_selection]
      val1 = val0[evs_selection, (order - 1):]
      val_fac = (np.arange(nEigvals-order+1, 0, -1).reshape((-1,1)) + np.arange(order-1).reshape((1, -1))).prod(axis=-1)
      val = np.concat([val1*val_fac, np.zeros((sum(evs_selection), order - 1), dtype=evs.dtype)], axis=1)
    ret[inds, :] = val
  return ret


def function_coeffs(M: np.ndarray, eigvals:np.ndarray|None=None):
  if eigvals is None:
    eigvals = np.linalg.eigvals(M)
  pass


t=np.arange(10)
M=np.zeros(t.shape*2)
M[t,t] = (1,2,2,3,3,4,4,4,4,4)
M[t[:-1],t[1:]] = (0,0,0,1,0,1,1,1,0)
# M = np.array(((2, 0), (0, 3)))
# M = np.array(((1, 1), (-1, 1)))
# M = np.array(((2, 1), (0, 2)))

M
evs = np.linalg.eigvals(M)
ev_mult = eigval_multiplicity(M, evs)
b = b_matrix(ev_mult)[:, ::-1].T
cs = scipy.linalg.solve(b, matrix_power_series(M, len(b)))


exp_M = np.tensordot(ev_mult.map(np.exp), cs, 1)
exp_M_ref = scipy.linalg.expm(M)
exp_M_ref2 = scipy.linalg.funm(M, np.exp)
sqrt_M = np.tensordot(ev_mult.map(np.sqrt), cs, 1)

np.linalg.norm(exp_M - exp_M_ref)
np.linalg.norm(exp_M - exp_M_ref2)
np.linalg.norm(sqrt_M @ sqrt_M - M)


log_M = np.tensordot(ev_mult.map(np.log), cs, 1)
np.linalg.norm(np.linalg.eigvals(log_M) - np.log(evs))


# See https://github.com/scipy/scipy/issues/21803#issuecomment-2455666759
big_M = np.array([[40, 1], [1, 0]])
big_evs = np.linalg.eigvalsh(big_M)
# big_evs = np.linalg.eigvals(big_M)
big_ev_mult = eigval_multiplicity(big_M, big_evs)
big_b = b_matrix(big_ev_mult)[:, ::-1].T
big_cs = scipy.linalg.solve(big_b, matrix_power_series(big_M, len(big_b)))
exp_big_M = np.tensordot(big_ev_mult.map(np.exp), big_cs, 1)
exp_big_M_ref = scipy.linalg.expm(big_M)
np.linalg.norm(exp_big_M - exp_big_M_ref)/np.sqrt(np.linalg.norm(exp_big_M)*np.linalg.norm(exp_big_M_ref))




