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
    ret = [fs[k](ev) for _, ev, _, k in self.ev_iter]
    return np.array(ret)

  @property
  def ev_iter(self):
    return ((i, ev, alg, k) for i, (ev, alg) in enumerate(zip(self.eigvals, self.algebraic)) for k in range(alg))


  @property
  def tr(self):
    return self.eigvals @ self.algebraic

  @property
  def det(self):
    return np.prod(self.eigvals ** self.algebraic)

  @property
  def dim(self):
    return np.sum(self.algebraic)

  @property
  def rank(self):
    return np.sum(self.algebraic - self.geometric + 1)

  @property
  def product_norm(self):
    if np.all(self.eigvals==0):
      return np.ones_like(self.eigvals[0])
    elif np.any(self.eigvals==0):
      inds = self.eigvals!=0
      return Multiplicity(*(q[inds] for q in self)).product_norm
    else:
      e = self.algebraic/np.sum(self.algebraic)
      return np.prod(np.abs(self.eigvals) ** e)


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


def function_coeffs(M: np.ndarray, eigvals:np.ndarray|None|Multiplicity=None, normalize_to:float|None=1.):
  # Idea behind the normalization: f(A) = f(sx) with x=A/s
  if isinstance(eigvals, Multiplicity):
    ev_mult = eigvals
  else:
    ev_mult = eigval_multiplicity(M, eigvals)
  if normalize_to is not None:
    norm_scale = ev_mult.product_norm*normalize_to**(1/ev_mult.dim)
    _ev_mult = Multiplicity(ev_mult.eigvals/norm_scale, *ev_mult[1:])
    _M = M/norm_scale
    cs_rescale = np.array([np.ones_like(norm_scale) if k==0 else norm_scale**k for _,_,_,k in ev_mult.ev_iter])
    b = b_matrix(_ev_mult)[:, ::-1].T
    _cs = scipy.linalg.solve(b, matrix_power_series(_M, len(b)))
    cs = cs_rescale.reshape((-1, 1, 1)) * _cs
  else:
    b = b_matrix(ev_mult)[:, ::-1].T
    cs = scipy.linalg.solve(b, matrix_power_series(M, len(b)))
  return cs, ev_mult

def apply_fn(M: np.ndarray, f, *dfs, gen_df=None, eigvals:np.ndarray|None|Multiplicity=None, normalize_to:float|None=1., **kwargs):
  if isinstance(f, str):
    f = f.lower().strip()
    if f == 'exp':
      f = np.exp
      gen_df = lambda _: np.exp
    elif f in 'log ln'.split():
      f = np.log
      gen_df = lambda k: (lambda x: np.prod(-np.arange(1, k))/x**k)
    elif f == 'inv':
      f = lambda x: x**-1
      gen_df = lambda k: (lambda x: np.prod(-np.arange(1, k+1))/x**(k+1))
    elif f == 'sin':
      f = np.sin
      _dfs = [np.sin, np.cos, lambda x: -np.sin(x), lambda x: -np.cos(x)]
      gen_df = lambda k: _dfs[k%4]
    elif f == 'cos':
      f = np.cos
      _dfs = [np.cos, lambda x: -np.sin(x), lambda x: -np.cos(x), np.sin]
      gen_df = lambda k: _dfs[k%4]
    elif f == 'sqrt':
      f = np.sqrt
      gen_df = lambda k: (lambda x: np.prod(np.arange(.5, 1-k, -1))/x**(k-.5))
    else:
      raise ValueError(f'Unknown function f={f}')
  cs, ev_mult = function_coeffs(M, eigvals, normalize_to)
  ret = np.tensordot(ev_mult.map(f, *dfs, gen_df=gen_df, **kwargs), cs, 1)
  return ret



t=np.arange(10)
M=np.zeros(t.shape*2)
M[t,t] = (1,2,2,3,3,4,4,4,4,4)
M[t[:-1],t[1:]] = (0,0,0,1,0,1,1,1,0)
# M = np.array(((2, 0), (0, 3)))
# M = np.array(((1, 1), (-1, 1)))
# M = np.array(((2, 1), (0, 2)))

M = np.random.randn(20, 20, 2) @ [1, 1j]

M
evs = np.linalg.eigvals(M)
ev_mult = eigval_multiplicity(M, evs)
b = b_matrix(ev_mult)[:, ::-1].T
cs = scipy.linalg.solve(b, matrix_power_series(M, len(b)))


exp_M = np.tensordot(ev_mult.map(np.exp), cs, 1)
exp_M_2 = np.tensordot(ev_mult.map(np.exp, gen_df=lambda _: np.exp), cs, 1)
exp_M_ref = scipy.linalg.expm(M)
sqrt_M = np.tensordot(ev_mult.map(np.sqrt), cs, 1)
sqrt_M_2 = apply_fn(M, 'sqrt')

def err(x,y):
  nxy = np.linalg.norm(x - y)
  nxy_s = np.sqrt(np.linalg.norm(x) * np.linalg.norm(y))
  return nxy if nxy_s==0 else nxy/nxy_s
  

# np.linalg.norm(exp_M - exp_M_ref)
# np.linalg.norm(exp_M_2 - exp_M_ref)
# np.linalg.norm(sqrt_M @ sqrt_M - M)
# np.linalg.norm(sqrt_M_2 @ sqrt_M_2 - M)
err(exp_M, exp_M_ref)
err(exp_M_2, exp_M_ref)
err(sqrt_M @ sqrt_M, M)
err(sqrt_M_2 @ sqrt_M_2, M)


log_M = apply_fn(M, 'log')
log_ev_mult = eigval_multiplicity(log_M)
log_ev_mult.product_norm
# np.linalg.norm(np.linalg.eigvals(log_M) - np.log(evs))
# np.linalg.norm(apply_fn(log_M, 'exp') - M)
# np.linalg.norm(apply_fn(2*log_M, 'exp') - M @ M)
# np.linalg.norm(apply_fn(10*log_M, 'exp') - matrix_power_series(M, 11)[-1])
err(np.linalg.eigvals(log_M), np.log(evs)) # Might be misleading for complex eigenvalues, also ordering
err(apply_fn(log_M, 'exp'), M)
err(apply_fn(2*log_M, 'exp'), M @ M)
err(apply_fn(10*log_M, 'exp'), matrix_power_series(M, 10 + 1)[-1])
err(apply_fn(100*log_M, 'exp'), matrix_power_series(M, 100 + 1)[-1])


# See https://github.com/scipy/scipy/issues/21803#issuecomment-2455666759
big_M = np.array([[40, 1], [1, 0]])
exp_big_M = apply_fn(big_M, 'exp')
exp_big_M_ref = scipy.linalg.expm(big_M)
err(exp_big_M, exp_big_M_ref)




