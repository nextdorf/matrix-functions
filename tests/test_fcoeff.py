from . import *
import scipy


frac_f_ref = lambda t: lambda A: scipy.linalg.fractional_matrix_power(A, t)

@pytest.mark.parametrize('args, real_vals, dim, count', direct_prod(
  [
    dict(f='exp', f_ref=scipy.linalg.expm),
    dict(f='log', f_ref=scipy.linalg.logm),
    dict(f='sin', f_ref=scipy.linalg.sinm),
    dict(f='sqrt', f_ref=scipy.linalg.sqrtm),
    dict(f='rnd_fractional', f_ref=frac_f_ref),
  ],
  (True, False),
  [2, 3, 4, 5, 7, 10, 15, 20, 30],
  5,
))
def test_apply_fn(args, real_vals: bool, dim: int, count: int):
  f, f_ref = map(args.get, 'f f_ref'.split())
  is_rnd_fractional = f == 'rnd_fractional'
  for _ in range(count):
    if is_rnd_fractional:
      t = np.random.randn()*5
      f = lambda x: (x + 0j)**t
      f_ref = frac_f_ref(t)
    if real_vals:
      M = np.random.randn(dim, dim)
    else:
      M = np.random.randn(dim, dim, 2) @ [1, 1j] / 2**.5
    f_M = apply_fn(M, f)
    f_ref_M = f_ref(M)
    same_nan = np.isnan(f_M) == np.isnan(f_ref_M)
    assert np.all(same_nan), repr(np.array([M, f_M, f_ref_M]))
    f_M = np.nan_to_num(f_M, nan=0)
    f_ref_M = np.nan_to_num(f_M, nan=0)
    assert err(f_M, f_ref_M) < 1e-9
