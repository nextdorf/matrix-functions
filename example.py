from matrixfuncs import *
import numpy as np

unit = .2
M = np.array([[2*np.cos(unit), -1], [1,0]])
krylovM = KrylovSpace(M)
krylovM.funcCoeffs.shape

mfn = MFunc.fromKrylov(krylovM)
mfn.coeffs
ts = np.linspace(0, 2*np.pi, 1000)
mfn0 = [1,0] @ mfn @ np.sin([0,-unit])
mfn0.coeffs
arr = mfn0(lambda t: np.exp(np.outer(np.log(t), ts/unit))).coeffs

import matplotlib.pyplot as plt
plt.plot(ts, arr)
plt.savefig('sinFromMFunc.png', dpi=200)



