import numpy as np
import scipy
from collections import namedtuple


def cluster(xs, zero_thrsh=1e-24, rel_thrsh=1e-8):
  xs = np.asarray(xs)
  assert len(xs.shape)==1, 'xs has to be a 1D array'
  zero_inds = abs(xs) < zero_thrsh
  has_zeros = np.any(zero_inds)
  non_zero_xs = xs[np.logical_not(zero_inds)]
  count_non_zero = len(non_zero_xs)
  if count_non_zero <= 1:
    pass

  non_zero_ks = np.arange(count_non_zero)
  for i in range(count_non_zero-1):
    for j in range(i+1, count_non_zero):
      rel_diff = abs(non_zero_xs[i]/non_zero_xs[j] - 1)
      if rel_diff < 2*rel_thrsh:
        non_zero_ks[j] = non_zero_ks[i]

  ks = set(non_zero_ks)
  ys = [np.mean(xs[k]) for k in ks]
  ys_count = [sum(non_zero_ks == k) for k in ks]
  if has_zeros:
    ys = [0] + ys
    ys_count = [sum(zero_inds)] + ys_count
  return np.asarray(ys), np.asarray(ys_count)


class Multiplicity(namedtuple('Multiplicity', 'eigvals algebraic geometric'.split())):
  def __new__(cls, eigvals, algebraic, geometric):
    eigvals = np.asarray(eigvals)
    algebraic = np.asarray(algebraic)
    geometric = np.asarray(geometric)
    assert len(eigvals.shape) == 1 and eigvals.shape == algebraic.shape == geometric.shape, 'Incompatible Shapes'
    return super(Multiplicity, cls).__new__(cls, eigvals, algebraic, geometric)

  @staticmethod
  def from_eigs(eigvals, eigvecs):
    eigvals = np.asarray(eigvals)
    eigvecs = np.asarray(eigvecs)
    evs, geom_mult = cluster(eigvals)
    ks = np.asarray([np.argmin(abs(evs - ev)) for ev in eigvals])
    alg_mult = []
    for k in range(len(evs)):
      eigspace = eigvecs[:, ks == k]
      span = scipy.linalg.orth(eigspace)
      alg_mult.append(span.shape[-1])
    return Multiplicity(evs, geom_mult, alg_mult)

  @staticmethod
  def from_matrix(M):
    return Multiplicity.from_eigs(*np.linalg.eig(M))


def funm(A, func, dfunc, eigval_multiplicity: Multiplicity = None):
  A: np.ndarray = np.asarray(A)
  if eigval_multiplicity is None:
    eigval_multiplicity = Multiplicity.from_matrix(A)
  evs = eigval_multiplicity.eigvals
  mult = eigval_multiplicity.algebraic - eigval_multiplicity.geometric + 1
  A_dim = A.shape[-1]
  dim = mult.sum()
  assert A.shape == (A_dim, A_dim), 'A has to be a square matrix'

  A0 = np.eye(A_dim, dtype=A.dtype, device=A.device)

  if dim == 1:
    return func(evs[0]) * A0

  at_non_zero_eigval = abs(evs) > 1e-16
  mult_at_non_zero_eigval = eigval_multiplicity.algebraic[at_non_zero_eigval]
  det_at_non_zero_eigval = (evs[at_non_zero_eigval]**mult_at_non_zero_eigval).prod()
  scale = abs(det_at_non_zero_eigval)**(1/mult_at_non_zero_eigval.sum())

  # By rescaling A such that all eigenvalues are of order 1 helps with stability
  A1 = A/scale
  normed_powers = [A0, A1]
  for i in range(2, dim):
    i_2 = i//2
    A_pow_i = normed_powers[i_2] @ normed_powers[i - i_2]
    normed_powers.append(A_pow_i)

  B = []
  func_vec = []
  for i, ev in enumerate(evs):
    B.append((ev/scale)**np.arange(dim))
    func_vec.append(func(ev))
    for j in range(1, mult[i]):
      b0 = B[-j][:-j]
      poch = (np.arange(1, dim-j+1).reshape((-1, 1)) + np.arange(j).reshape((1, -1))).prod(axis=1)
      b1 = np.concat([(0,)*j, b0*poch])
      B.append(b1)
      func_vec.append(dfunc(j, ev) * scale**j)
  B = np.asarray(B).swapaxes(0, 1)

  coeffs = scipy.linalg.solve(B, normed_powers)
  func_A = np.tensordot(func_vec, coeffs, 1)
  return func_A

