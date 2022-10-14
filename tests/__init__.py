from matrixfuncs import *
from matrixfuncs.utils import *
import numpy as np
import pytest

def directProd(xs0, *xss):
  if len(xss) == 0:
    for x in xs0:
      yield (x,)
  else:
    inner = list(directProd(*xss))
    for x in xs0:
      for y in inner:
        yield (x,) + y

def rndIndexPermutation(n: int):
  return IndexPermutation(*np.random.choice(n, (n,), False))

def rndPermutation(m: int, n: int):
  vals = np.random.choice(m, (n,))
  return Permutation(vals, normalOrder=sorted(vals))


