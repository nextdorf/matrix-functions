from . import *
import sympy as sp


def assertSympyEquality(a, b, msg = ''):
  assert (np.array((sp.Array(a) - sp.Array(b)).simplify().tolist()) == 0).all(), msg

def assertSymKrylovSpaceEquality(sk1, sk2):
  for attrName in 'lambdaVec betaVec betaM lambdaCoeffs funcCoeffs'.split():
    sk1Attr = getattr(sk1, attrName)
    sk2Attr = getattr(sk2, attrName)
    assertSympyEquality(sk1Attr, sk2Attr, f'{attrName} differ')
  assertSympyEquality(sk1.eigvals, sk2.eigvals, 'eigvals differ')
  assert (sk1.multiplicity-sk2.multiplicity == 0).all(), 'multiplicity differ'

class TestSymKrylovSpace:
  cache = SymKrylovSpaceCache()

  @pytest.mark.parametrize('multiplicity', [[1], [1,2], [2,1], [1,3,2], [3,2,1], [1,2,3]])
  def test_cache(self, multiplicity):
    sk = SymKrylovSpace(multiplicity)
    csk = TestSymKrylovSpace.cache[multiplicity]
    assertSymKrylovSpaceEquality(sk, csk)

  @pytest.mark.parametrize('multiplicity', [[1], [1,2], [2,1], [3,2], [4,1]])#, [2,1,1]])
  def test_dimensionReduction(self, multiplicity: list[int]):
    sk = TestSymKrylovSpace.cache[multiplicity]
    for i,_ in enumerate(multiplicity):
      fcRed = sk.reduceFunctionCoeffs(*np.array(sk.eigvals)[[i, i]])
      assertSympyEquality(sk.funcCoeffs, fcRed, repr((i, i)))
      for j,_ in list(enumerate(multiplicity))[i+1:]:
        multIJ = multiplicity[:i]+[multiplicity[i]+multiplicity[j]]+multiplicity[i+1:j]+multiplicity[j+1:]
        fcij = TestSymKrylovSpace.cache[multIJ].funcCoeffs
        fcRedij = sk.reduceFunctionCoeffs(*np.array(sk.eigvals)[[i, j]])
        assertSympyEquality(fcij, fcRedij, repr((i, j)))
        fcRedji = sk.reduceFunctionCoeffs(*np.array(sk.eigvals)[[j, i]])
        assertSympyEquality(fcij, fcRedji, repr((j, i)))

    #assertSymKrylovSpaceEquality(sk, csk)


