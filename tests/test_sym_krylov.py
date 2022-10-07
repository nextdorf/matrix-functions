from matrixfuncs import *
import sympy as sp
import numpy as np

def test_case0():
  c = SymKrylovSpaceCache()
  multiplicity = [1,2]

  sk = SymKrylovSpace(multiplicity)
  csk = c[multiplicity]

  for attrName in 'eigvals lambdaVec betaVec betaM lambdaCoeffs funcCoeffs'.split():
    attr: sp.Array = getattr(sk, attrName)
    attrC: sp.Array = getattr(csk, attrName)
    isZero = np.all(np.array((attr - attrC).simplify().tolist()) == 0)
    assert isZero, f'{attrName}: Cached evaluation and direct evaluation differ'
