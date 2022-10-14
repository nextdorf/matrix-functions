from typing import Generator
from . import *


class TestIndexPermutation:
  @pytest.mark.parametrize('n', range(11))
  def test_One(self, n: int):
    permOne = IndexPermutation.One(n)
    perm = np.array(permOne)
    assert np.all(perm == np.arange(n))
    for _ in range(n):
      permRnd = rndIndexPermutation(n)
      perm = np.array(permRnd)
      perm1 = np.array(permRnd * permOne)
      perm2 = np.array(permOne * permRnd)
      assert np.all(perm1 == perm)
      assert np.all(perm2 == perm)

  @pytest.mark.parametrize('perm', map(rndIndexPermutation, [1,2,3,10,20, 50, 100, 1000, 10000]))
  def test_Inverse(self, perm: IndexPermutation):
    perm1 = np.array(perm * perm.inverse)
    perm2 = np.array(IndexPermutation.One(len(perm)))
    assert np.all(perm1 == perm2)

  @pytest.mark.parametrize('perm, n', directProd(
      map(rndIndexPermutation, [1,2,3,10,20,30,50,70,80,100,200]),
      [0,-1, None, 2, 10, 33, 100]
    ))
  def test_Pow(self, perm: IndexPermutation, n: int|None):
    if n is None:
      n = len(perm)
    perm1 = np.array(perm ** n)
    if n >= 0:
      perm2 = IndexPermutation.One(len(perm))
      for _ in range(n):
        perm2 *= perm
    else:
      perm2 = IndexPermutation.One(len(perm))
      permInv = perm.inverse
      for _ in range(-n):
        perm2 *= permInv
    assert np.all(perm1 == np.array(perm2))

  @pytest.mark.parametrize('x, y', map(
      lambda q: (q[0](np.random.choice(q[1], (len(q[2]),))), q[2]),
      directProd(
        [list, lambda xs: (x for x in xs), tuple],
        [1,2,3,10,20,50,100],
        map(rndIndexPermutation, [1,2,3,10,20,50,100])
    )))
  def test_matmul(self, x: list|Generator|tuple, y: IndexPermutation):
    x_y = x @ y
    if any([isinstance(x, t) for t in [list, Generator, tuple]]):
      assert isinstance(x_y, type(x))
    else:
      assert isinstance(x_y, list)
    if not isinstance(x, Generator): #x is exausted otherwise
      assert np.all(np.array(list(x_y)) == np.array(list(x))[list(y)])

class TestPermutation:
  @pytest.mark.parametrize('n', [1,2,3,10,20, 50, 100, 1000])
  def test_findIndexPermutation(self, n: int):
    vals = list(np.random.choice(10*n, (n,)))
    perm0 = rndIndexPermutation(len(vals))
    normalOrder: list = vals @ perm0.inverse
    perm1 = Permutation.findIndexPermutation(vals, normalOrder)
    vals1 = normalOrder @ perm1
    assert np.all(np.array(vals) == np.array(vals1))



