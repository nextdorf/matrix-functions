from . import *
import sympy as sp

x, y, z = sp.var('x y z')
f = sp.Function('f')

@pytest.mark.parametrize('expr, res, n', [(
    f(y)*(x/(x-y))**2 + x*y*f(x).diff()/(x-y) - y*(2*x-y)*f(x)/(x-y)**2,
    f(x) - x*f(x).diff() + x**2*f(x).diff(x, 2)/2,
    2 ),
  ])
def test_fn(expr: sp.Expr, res: sp.Expr, n:int):
  assert np.all(np.array((differentialLimit(expr, y, x, n)-res).simplify()) == 0), 'y -> x'
  assert np.all(np.array((differentialLimit(expr, x, y, n)-res.subs(x,y)).simplify()) == 0), 'x -> y'


