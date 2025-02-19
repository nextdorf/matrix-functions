import numpy as np
from collections import namedtuple
import numdifftools as nd
import scipy

class Multiplicity(namedtuple('Multiplicity', 'eigvals algebraic geometric'.split())):
  '''A named tuple representing the multiplicities of eigenvalues in a matrix.

  This class stores distinct eigenvalues along with their algebraic and geometric multiplicities.
  It provides methods for applying functions to eigenvalues and computing properties like trace,
  determinant, rank, and a normalization factor.

  Attributes
  ----------
  eigvals : np.ndarray
    An array containing distinct eigenvalues of the matrix.
  algebraic : np.ndarray
    The algebraic multiplicity of each corresponding eigenvalue.
  geometric : np.ndarray
    The geometric multiplicity of each corresponding eigenvalue.
  '''
  def __new__(cls, eigvals: np.ndarray, algebraic: np.ndarray, geometric: np.ndarray):
    '''Constructor for the Multiplicity class.

    Parameters
    ----------
    eigvals : np.ndarray
      Distinct eigenvalues of the matrix.
    algebraic : np.ndarray
      Algebraic multiplicities of each eigenvalue.
    geometric : np.ndarray
      Geometric multiplicities of each eigenvalue.

    Returns
    -------
    Multiplicity
      A named tuple containing eigenvalues and their multiplicities.
    '''
    return super(Multiplicity, cls).__new__(cls, np.asarray(eigvals), np.asarray(algebraic), np.asarray(geometric))
  
  def map(self, f, *dfs, gen_df=None, mult_space='min', **kwargs):
    '''Maps a function (and derivatives) to the eigenvalues.

    This method applies `f` and its derivatives to each eigenvalue according to its
    algebraic multiplicity.

    Parameters
    ----------
    f : function
      The function to apply to the eigenvalues.
    dfs : tuple of functions, optional
      Precomputed derivatives of `f`, if available.
    gen_df : function, optional
      A function that generates higher-order derivatives when needed. gen_df(k) should generate
      the kth derivative. Using dfs takes precedance over using gen_df. If gen_df is `None` then
      `numdifftools.Derivative` is used.

    Returns
    -------
    ndarray
      (f(ev0), f'(ev0), f''(ev0), ..., f(ev1), f'(ev1), f''(ev1), ...)
      
    Notes
    -----
    See also `ev_iter` for further details on the order.
    '''
    diff_count = np.max(self.algebraic)
    fs = [f] + list(dfs[:diff_count - 1])

    if len(fs) < diff_count:
      if gen_df is None:
        k0 = len(fs) - 1
        gen_df = lambda i, **_: nd.Derivative(fs[k0], n=i - k0, order=2 * 1)
      for i in range(len(fs), diff_count):
        fs.append(gen_df(i))

    #TODO: vectorize
    ev_iter = self.ev_iter(mult_space)
    ret = [fs[k](ev) for _, ev, _, k in ev_iter]
    return np.array(ret)

  def ev_iter(self, mult_space='min'):
    '''Iterate over eigenvalues and their multiplicities.

    Parameters
    ----------
    mult_space : str
      TODO

    Returns
    -------
    generator
      Yields (eigenvalue index, eigenvalue, multiplicity, counter of the multiplicity)

      The inner iteration is through the multiplicity and the outer iterator is through the
      eigenvalues.
    '''
    return ((i, ev, mult, k) for i, (ev, mult) in enumerate(zip(self.eigvals, self._multiplicity(mult_space))) for k in range(mult))

  def _multiplicity(self, mult_space='min'):
    if mult_space == 'min':
      return self.algebraic - self.geometric + 1
    elif mult_space == 'full':
      return self.algebraic
    else:
      raise ValueError(f'Unsupported mult_space "{mult_space}"')


  @property
  def tr(self):
    'The trace'
    return self.eigvals @ self.algebraic

  @property
  def full_eigvals(self):
    'All eigenvalues'
    return np.array([ev for _,ev,_,_ in self.ev_iter('full')])

  @property
  def det(self):
    'The determinant'
    return np.prod(self.eigvals ** self.algebraic)

  @property
  def dim(self):
    'The dimension'
    return np.sum(self.algebraic)

  @property
  def rank(self):
    'The rank'
    return np.sum(self.algebraic - self.geometric + 1)

  def product_norm(self, mult_space='min'):
    '''Compute a norm-like quantity used for matrix normalization.

    - If all eigenvalues are nonzero, returns `|det|^(1/dim)`.
    - If some but not all eigenvalues are nonzero, returns `product_norm` of the nonzero part.
    - If the only eigenvalue is 0, returns 1.

    Returns
    -------
    float
      Computed normalization value.
    '''
    if np.all(self.eigvals == 0):
      return np.ones_like(self.eigvals[0])
    elif np.any(self.eigvals == 0):
      inds = self.eigvals != 0
      return Multiplicity(*(q[inds] for q in self)).product_norm(mult_space)
    else:
      mult = self._multiplicity(mult_space)
      e = mult / np.sum(mult)
      return np.prod(np.abs(self.eigvals) ** e)

  @staticmethod
  def from_matrix(M: np.ndarray, eigvals: np.ndarray | None = None, **kwargs):
    '''Create a `Multiplicity` instance from a given matrix.

    Parameters
    ----------
    M : np.ndarray
      The input matrix.
    eigvals : np.ndarray, optional
      Precomputed eigenvalues. If `None`, they are computed automatically.

    Returns
    -------
    Multiplicity
      A `Multiplicity` instance containing eigenvalue data.

    See also
    --------
    `eigval_multiplicity`
    '''
    return eigval_multiplicity(M=M, eigvals=eigvals, **kwargs)


def matrix_power_series(M: np.ndarray, stop: int):
  '''Compute a sequence of matrix powers up to a given order.

  This function efficiently computes powers of a square matrix up to `stop - 1` using recursive
  squaring to reduce computational complexity.

  Parameters
  ----------
  M : np.ndarray
    A square matrix whose powers need to be computed.
  stop : int
    The exclusive upper bound to compute in the series.

  Returns
  -------
  np.ndarray
    An array containing matrices from the identity matrix (I) to M^(stop-1).
  '''
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
  '''Compute the algebraic and geometric multiplicities of eigenvalues.

  This function determines the distinct eigenvalues of a given matrix, their algebraic
  multiplicities, and geometric multiplicities (dimension of the eigenspaces).

  Parameters
  ----------
  M : np.ndarray
    A square matrix whose eigenvalues are analyzed.
  eigvals : np.ndarray, optional
    Precomputed eigenvalues of the matrix. If None, they are computed internally.
  zero_thrsh : float = 1e-15
    Threshold below which eigenvalues are considered zero. Set to 0 if only 0 itself should be
    considered a vanishing eigenvalue.
  rel_eq_thrsh : float = 1e-8
    Relative threshold for treating eigenvalues as identical.

  Returns
  -------
  mult: Multiplicity
  '''
  # TODO: This seems to be extremely instable for eigenvalues which are close but not equal, however fcoeffs.dim seems to do a good job in finding the dimension
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

def b_matrix(multiplicity: Multiplicity, mult_space='min'):
  '''Construct the basis matrix `phi` for function approximation using eigenvalue decomposition.

  This function generates a matrix used to solve for function coefficients when applying matrix
  functions.

  Parameters
  ----------
  multiplicity : Multiplicity
    Eigenvalue multiplicity structure computed from a matrix.

  Returns
  -------
  np.ndarray
    A transformation matrix for function coefficient computation. Intended to be used like
    `cs = scipy.linalg.solve(b, matrix_power_series(M, len(b)))`. See also the example.

  Example
  -------
  >>> # Construct a random 5x5 matrix M and explicitly construct the exp(M)
  >>> M = np.random.randn(5,5)
  >>> mult = eigval_multiplicity(M)
  >>> b = b_matrix(mult)[:, ::-1].T
  >>> cs = scipy.linalg.solve(b, matrix_power_series(M, len(b)))
  >>> 
  >>> # M can be set to None since `coeffs` and `eigvals` are provided
  >>> expM = apply_fn(None, 'exp', coeffs=cs, eigvals=mult)
  >>> 
  >>> # Compare solution with reference solution from `scipy.linalg.expm`
  >>> expM_ref = scipy.linalg.expm(M)
  >>> np.linalg.norm(expM - expM_ref)
  >>>
  >>> #The expected rounding error is of order 1.11e-16 * order of biggest values * sqrt(number of calculation steps in that biggest order)
  np.float64(2.0755569795398968e-15)
  '''
  ret = []
  ev_iter = list(multiplicity.ev_iter(mult_space))
  dim = len(ev_iter)
  for (_, ev, mult, k) in ev_iter:
    val0 = np.concat([ev**np.arange(dim-1-k, 0, -1), [1]])
    val_fac = (np.arange(dim-k, 0, -1).reshape((-1, 1)) + np.arange(k).reshape((1, -1))).prod(axis=-1)
    val = np.concat([val0 * val_fac, [0]*k])
    ret.append(val)
  ret = np.array(ret)
  return ret

def function_coeffs(M: np.ndarray, eigvals:np.ndarray|None|Multiplicity=None, mult_space='min', normalize_to:float|None=1.):
  '''Compute coefficients for matrix function computation.

  This function determines the necessary coefficients to compute functions applied to matrices,
  leveraging eigenvalue decomposition.

  Parameters
  ----------
  M : np.ndarray
    The input square matrix.
  eigvals : np.ndarray or Multiplicity, optional
    Precomputed eigenvalues or their multiplicities. If it is not of type Multiplicity then it is
    calculated using `eigval_multiplicity`
  normalize_to : float, optional, default=1.
    Scaling factor to normalize eigenvalues, improving numerical stability.

    To achieve better accuracy the `M` is rescaled such that the eigenvalues are approximately of
    size `normalize_to`. In particular, for a non-singular Matrix its determinate will be rescaled
    to `normalize_to`. Setting `normalize_to` to `None` skips the normalization step.

  Returns
  -------
  cs,mult : ndarray, Multiplicity
    The quantities necessary for applying a function to a matrix. `cs` is the `phi`-tensor. `mult`
    is the the multiplicity used. If the type of `eigvals` is `Multiplicity` then `mult` is equal
    to `eigvals`.

  See also
  --------
  See also `b_matrix` or `apply_fn` for an example
  '''

  # Idea behind the normalization: f(A) = f(sx) with x=A/s
  #TODO: Allow for reduced multiplicity
  if isinstance(eigvals, Multiplicity):
    ev_mult = eigvals
  else:
    ev_mult = eigval_multiplicity(M, eigvals)
  if normalize_to is not None:
    ev_iter = list(ev_mult.ev_iter(mult_space))
    norm_scale = ev_mult.product_norm(mult_space)/normalize_to**(1/len(ev_iter))
    _ev_mult = Multiplicity(ev_mult.eigvals/norm_scale, *ev_mult[1:])
    _M = M/norm_scale
    cs_rescale = np.array([np.ones_like(norm_scale) if k==0 else norm_scale**k for _,_,_,k in ev_iter])
    b = b_matrix(_ev_mult, mult_space=mult_space)[:, ::-1].T
    _cs = scipy.linalg.solve(b, matrix_power_series(_M, len(b)))
    cs = cs_rescale.reshape((-1, 1, 1)) * _cs
  else:
    b = b_matrix(ev_mult, mult_space=mult_space)[:, ::-1].T
    cs = scipy.linalg.solve(b, matrix_power_series(M, len(b)))
  return cs, ev_mult

def apply_fn(M: np.ndarray, f, *dfs, gen_df=None, eigvals:np.ndarray|None|Multiplicity=None, coeffs:np.ndarray|None=None, normalize_to:float|None=1., **kwargs):
  '''Apply a scalar function to a matrix using spectral decomposition.

  If `eigvals` is of type Multiplicity and `coeffs` is provided then `M` is ignored and the
  function is applied directly.

  Parameters
  ----------
  M : np.ndarray
    The square matrix to transform.
  f : function
    A function which is applied on the eigenvalues. Can be a predefined function like 'exp', 'log',
    etc. If a supported function is provided `dfs` and `gen_df` are ignored.

    Full list of supported predefined functions:
    - `exp`
    - `log` or `ln`
    - `inv`
    - `sin`, and `cos`
    - `sqrt`
  dfs : function
    Precomputed derivatives of `f`, if available.
  gen_df : function, optional
    A function that generates higher-order derivatives when needed. gen_df(k) should generate
    the kth derivative. Using dfs takes precedance over using gen_df. If gen_df is `None` then
    `numdifftools.Derivative` is used.
  eigvals : Optional[ndarray|Multiplicity] = None
    Precomputed eigenvalues or their multiplicities. If it is not of type Multiplicity then it is
    calculated using `eigval_multiplicity`
  cs : Optional[ndarray]
    The `phi`-tensor. See also `function_coeffs`
  normalize_to : Optional[float] = 1.
    Scaling factor to normalize eigenvalues, improving numerical stability.

    To achieve better accuracy the `M` is rescaled such that the eigenvalues are approximately of
    size `normalize_to`. In particular, for a non-singular Matrix its determinate will be rescaled
    to `normalize_to`. Setting `normalize_to` to `None` skips the normalization step.

  Returns
  -------
  f_M: ndarray
    The matrix after applying the function `f`.

  Example
  -------
  >>> # Construct a random 5x5 matrix M and explicitly construct the exp(M)
  >>> # Compare solution with reference solution from `scipy.linalg.expm`
  >>> M = np.random.randn(5,5)
  >>> expM = apply_fn(M, 'exp')
  >>> expM_ref = scipy.linalg.expm(M)
  >>> np.linalg.norm(expM - expM_ref)
  >>>
  >>> #The expected rounding error is of order 1.11e-16 * order of biggest values * sqrt(number of calculation steps in that biggest order)
  np.float64(2.0755569795398968e-15)

    See also
  --------
  See also `b_matrix` for more fine-grained usage, which might be preferred if many functions
  have to be calculated for the same matrix.
  '''
  if isinstance(f, str):
    f = f.lower().strip()
    if f == 'exp':
      f = np.exp
      gen_df = lambda _: np.exp
    elif f in 'log ln'.split():
      _0 = np.array(0j)
      f = lambda x: np.log(x + _0)
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
      _0 = np.array(0j)
      f = lambda x: np.sqrt(x + _0)
      gen_df = lambda k: (lambda x: np.prod(np.arange(.5, 1-k, -1))/x**(k-.5))
    else:
      raise ValueError(f'Unknown function f={f}')
  if coeffs is not None and isinstance(eigvals, Multiplicity):
    cs, ev_mult = coeffs, eigvals
  else:
    _kwargs = {k: kwargs[k] for k in 'mult_space'.split() if k in kwargs}
    cs, ev_mult = function_coeffs(M, eigvals, normalize_to=normalize_to, **_kwargs)
  f_ev_mult = ev_mult.map(f, *dfs, gen_df=gen_df, **kwargs)
  ret = np.tensordot(f_ev_mult, cs, ((0, ), (0, )))
  return ret
