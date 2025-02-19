import numpy as np

def err(x,y):
  nxy = np.linalg.norm(x - y)
  nxy_s = np.sqrt(np.linalg.norm(x) * np.linalg.norm(y))
  return nxy if nxy_s==0 else nxy/nxy_s
  

