import numpy as np

def err(x,y):
  '''Compute the relative error between two matrices or vectors.

  Parameters
  ----------
  x : np.ndarray
    First input array (matrix or vector).
  y : np.ndarray
    Second input array (matrix or vector) to compare against `x`.
  
  Returns
  -------
  float
      The relative error metric, or the absolute norm difference if both `x` and `y` have zero norm.
  '''
  nxy = np.linalg.norm(x - y)
  nxy_s = np.sqrt(np.linalg.norm(x) * np.linalg.norm(y))
  return nxy if nxy_s==0 else nxy/nxy_s
  

