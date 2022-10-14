from typing import Generator
import numpy as np

class IndexPermutation:
  def __init__(self, *values: int, _assertCorrectness=True):
    self.value = tuple(values)
    if _assertCorrectness:
      assert np.all(np.sort(self.value) == np.arange(len(self.value))), 'Is not a permutation'

  @property
  def inverse(self):
    _, invValue = zip(*sorted(zip(self.value, range(len(self))), key=lambda x: x[0]))
    return IndexPermutation(*invValue, _assertCorrectness=False)

  @staticmethod
  def One(n: int):
    return IndexPermutation(*range(n), _assertCorrectness=False)

  def swap(self, i: int, j: int):
    n = len(self)
    if i%n==j%n: return self
    val = list(self.value)
    vi = val[i]
    val[i] = val[j]
    val[j] = vi
    return IndexPermutation(*val, _assertCorrectness=False)

  def __mul__(self, other: 'IndexPermutation'):
    if isinstance(other, IndexPermutation):
      pass
    elif hasattr(other, '__iter__'):
      other = IndexPermutation(*other)
    else:
      return NotImplemented
    if len(self) != len(other):
      return NotImplemented
    return IndexPermutation(*(self[i] for i in other), _assertCorrectness=False)

  def __matmul__(self, other):
    if isinstance(other, IndexPermutation):
      return self * other
    else:
      return NotImplemented

  def __rmatmul__(self, other):
    otherT = type(other)
    if otherT is IndexPermutation:
      return self * other
    elif not hasattr(other, '__iter__'):
      return NotImplemented
    x = list(other)
    ret = (x[i] for i in self)
    if  isinstance(other, tuple): return tuple(ret)
    elif  isinstance(other, Generator): return ret
    else: return list(ret)

  def __pow__(self, other: int):
    if not isinstance(other, int):
      return NotImplemented
    '''
    n = len(self)-1
    m = other % n
    if m > 0 and m - n/2 > 0:
      m -= n

    if m==0:
      return IndexPermutation.One(n+1)
    elif m>0:
      perm, m = self, m
    else:
      perm, m = self.inverse, -m

    ret = perm
    while m > 1:
      ret = ret*perm
      m-=1
    return ret
    '''
    if other==0:
      return IndexPermutation.One(len(self))
    elif other>0:
      perm, m = self, other
    else:
      perm, m = self.inverse, -other

    ret = perm
    while m > 1:
      ret = ret*perm
      m-=1
    return ret

  def __getitem__(self, idx):
    return self.value[idx]
  def __len__(self):
    return len(self.value)
  def __iter__(self):
    return iter(self.value)

  def __str__(self) -> str:
    return str(self.value)
  def __repr__(self) -> str:
    s0 = str(self)[1 : -1]
    return f'IndexPermutation[{s0}]'


class Permutation:
  def __init__(self, value: list, perm: IndexPermutation|None=None, normalOrder: list|None = None, _assertCorrectness = True):
    self.value = list(value)
    if perm is not None:
      self.perm = perm
      if normalOrder is not None:
        self.normalOrder = list(normalOrder)
      else:
        self.normalOrder = self.value @ self.perm.inverse
    else:
      if normalOrder is not None:
        self.normalOrder = normalOrder
        self.perm = Permutation.findIndexPermutation(self.value, self.normalOrder, _assertCorrectness=_assertCorrectness)
      else:
        self.perm = IndexPermutation.One(len(self.value))
        self.normalOrder = self.value

  @property
  def inverse(self):
    inv = self.perm.inverse
    return Permutation(self.normalOrder @ inv, inv, self.normalOrder, False)

  @staticmethod
  def findIndexPermutation(value: list, normalOrder: list, _assertCorrectness = True):
    if _assertCorrectness:
      assert sorted(value) == sorted(normalOrder), 'value is not a permutation of normalOrder'
    permOfValue = []
    if not hasattr(value, '__len__'):
      value = list(value)
    if not hasattr(normalOrder, '__getitem__'):
      normalOrder = list(normalOrder)
    rs = list(range(len(value)))
    for v in value:
      for idx, i in enumerate(rs):
        if normalOrder[i] == v:
          permOfValue.append(i)
          del rs[idx]
          break
    return IndexPermutation(*permOfValue, _assertCorrectness=_assertCorrectness)


  def swap(self, i: int, j: int):
    if i==j:
      return self**1
    else:
      return self @ IndexPermutation.One(len(self)).swap(i, j)

  def __matmul__(self, other: 'IndexPermutation | Permutation'):
    if isinstance(other, IndexPermutation):
      return Permutation(self.value @ other, self.perm @ other, self.normalOrder, False)
    elif isinstance(other, Permutation) and self.normalOrder==other.normalOrder:
      return self @ other.perm
    else:
      return NotImplemented

  def __pow__(self, other: int):
    if not isinstance(other, int):
      return NotImplemented
    if other == 0:
      return Permutation(self.normalOrder, IndexPermutation.One(len(self)), self.normalOrder, False)
    elif other == 1:
      return Permutation(self.value, self.perm, self.normalOrder, False)
    elif other == -1:
      return self.inverse
    perm = self.perm**other
    return Permutation(self.normalOrder @ perm, perm, self.normalOrder, False)


  def __getitem__(self, idx):
    return self.value[idx]
  def __len__(self):
    return len(self.value)
  def __iter__(self):
    return iter(self.value)

  def __str__(self) -> str:
    return f'Permutation[{str(self.value)[1:-1]} | {str(self.perm.value)[1:-1]}]' if len(self) else 'Permutation[]'
  def __repr__(self) -> str:
    return str(self)

