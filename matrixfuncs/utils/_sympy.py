import sympy as sp

def differentialLimit(expr: sp.Expr, x: sp.Symbol, y: sp.Symbol, n=2):
  '''Takes the limit of `expr` at `x` -> `y` but ensures that difference quotients are correctly replaced by derivatives.
  '''
  def parseName(gnName):
    idStr = '^{D:'
    idIdx = gnName.rfind(idStr)
    if idIdx < 0: return None, None
    fnName = gnName[:idIdx]
    idIdx1 = gnName.find('}', idIdx)
    fnArgs = gnName[idIdx+len(idStr):idIdx1][1:-1].strip()
    fnDiffs = []
    i0 = fnArgs.find('(', 0)
    while i0 >= 0:
      i1 = fnArgs.find(')', i0)
      sVar, sInt = map(lambda s: s.strip(), fnArgs[1:i1].split(','))
      xVar, xInt = sp.var(sVar), int(sInt)
      if xVar == x: xVar = y
      fnDiffs.append((xVar, xInt))
      i0 = fnArgs.find('(', i1)
    return fnName, fnDiffs
  def getDiffMapping(d: sp.Derivative):
    (fn, *args) = d.args
    fnName = str(fn.func)
    fnArgs = fn.args
    gn = sp.Function(f'{fnName}^{{D:{list(args)}}}')(*fnArgs)
    return gn
  def invDiffMapping(d: sp.Derivative | sp.Function):
    args, gnName, gnArgs = [None]*3
    if isinstance(d, sp.Derivative):
      (gn, *args) = d.args
      gnName = str(gn.func)
      gnArgs = gn.args
    elif isinstance(d, sp.Function):
      args = []
      gnName = str(d.func)
      gnArgs = d.args
    fnName, fnDiffs = parseName(gnName)
    if fnName is None: return None
    fn: sp.Derivative|sp.Function = sp.Function(fnName)(*gnArgs).diff(*fnDiffs, *args)
    return fn
  eps = sp.Dummy('\\epsilon')
  diffs = expr.atoms(sp.Derivative)
  diffMap = {d: dm for d,dm in ((_d, getDiffMapping(_d)) for _d in diffs) if dm is not None}
  expr0: sp.Expr = sp.series(expr.subs(list(diffMap.items())).subs(x,y+eps), eps, 0, n).doit()

  diffs1 = expr0.atoms(sp.Derivative, sp.Function)
  diffMap1 = {d: dm for d,dm in ((_d, invDiffMapping(_d)) for _d in diffs1) if dm is not None}
  expr1 = expr0.subs(list(diffMap1.items()))
  expr2 = sp.limit(expr1, eps, 0).doit()
  return expr2

