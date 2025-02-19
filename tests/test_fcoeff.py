from . import *
import scipy


@pytest.mark.parametrize('args, real_vals, dim, count', direct_prod(
  [dict(f=q[0], f_ref=q[1], dfs=(), gen_df=None) for q in (
    ('exp', scipy.linalg.expm),
    ('log', scipy.linalg.logm),
    ('sin', scipy.linalg.sinm),
    ('sqrt', scipy.linalg.sqrtm),
  )],
  (True, False),
  [2, 3, 4, 5, 7, 10, 15, 20, 30],
  5,
))
def test_gl(args, real_vals: bool, dim: int, count: int):
  f, dfs, gen_df, f_ref = map(args.get, 'f dfs gen_df f_ref'.split())
  for _ in range(count):
    if real_vals:
      M = np.random.randn(dim, dim)
    else:
      M = np.random.randn(dim, dim, 2) @ [1, 1j] / 2**.5
    f_M = apply_fn(M, f, *dfs, gen_df=gen_df)
    f_ref_M = f_ref(M)
    same_nan = np.isnan(f_M) == np.isnan(f_ref_M)
    assert np.all(same_nan), repr(np.array([M, f_M, f_ref_M]))
    f_M = np.nan_to_num(f_M, nan=0)
    f_ref_M = np.nan_to_num(f_M, nan=0)
    assert err(f_M, f_ref_M) < 1e-9
