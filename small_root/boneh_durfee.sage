from sage.all import *
from copy import copy
import itertools
import sys

# display matrix picture with 0 and X
# references: https://github.com/mimoo/RSA-and-LLL-attacks/blob/master/boneh_durfee.sage
def matrix_overview(BB):
  for ii in range(BB.dimensions()[0]):
    a = ('%02d ' % ii)
    for jj in range(BB.dimensions()[1]):
      a += ' ' if BB[ii,jj] == 0 else 'X'
      if BB.dimensions()[0] < 60:
        a += ' '
    print a
  print 

def boneh_durfee_bivariate_wiener_bound(_pol, modulo, XX, YY, mm, tt):
  '''
  Boneh-Durfee's Algorithm for Solving Modular Bivariate Equation with Small Root - Wiener's Bound
  '''
  PR = PolynomialRing(ZZ, 'x, y', order='degneglex')
  x, y = PR.gens()
  pol = PR(_pol)
  gik = {}
  hjk = {}
  for k in xrange(mm + 1):
    fk = pol^k * modulo^(mm-k)
    for i in xrange(mm - k + 1):
      gik[i, k] = x^i * fk
    for j in xrange(1, tt + 1):
      hjk[j, k] = y^j * fk

  monomials = set()
  for k in gik.keys():
    monomials |= set(gik[k].monomials())
  for k in hjk.keys():
    monomials |= set(hjk[k].monomials())
  monomials = list(monomials)
  monomials.sort()
  M = Matrix(ZZ, ((mm + 1) * (mm + 2)) / 2 + (mm + 1) * tt)

  row = 0
  for k in xrange(mm + 1):
    for i in xrange(mm - k + 1):
      for m, monomial in enumerate(monomials):
        M[row, m] = gik[i, k].monomial_coefficient(monomial) * monomial(XX, YY)
      row += 1
  for k in xrange(mm + 1):
    for j in xrange(1, tt + 1):
      for m, monomial in enumerate(monomials):
        M[row, m] = hjk[j, k].monomial_coefficient(monomial) * monomial(XX, YY)
      row += 1

  matrix_overview(M)

  B = M.LLL()

  matrix_overview(B)

  Hi = {}
  for i in xrange(B.nrows()):
    Hi[i] = 0
    for m, monomial in enumerate(monomials):
      Hi[i] += ZZ(B[i, m] / monomial(XX, YY)) * monomial

  PX = PolynomialRing(ZZ, 'xk')
  PY = PolynomialRing(ZZ, 'yk')

  xk = PX.gen()
  yk = PY.gen()

  root_x = root_y = None

  for i in xrange(len(Hi)):
    for j in xrange(len(Hi)):
      if i == j:
        continue
      pol1 = Hi[i]
      pol2 = Hi[j]
      pol_res = pol1.resultant(pol2, y).subs(x=xk)
      if not isinstance(pol_res, Integer):
        roots = pol_res.roots()
        if len(roots) > 0:
          if len(roots) == 1 and roots[0][0] == 0:
            continue
          _root_x = filter(lambda t: 0 < abs(t[0]) < XX, roots)
          if len(_root_x) == 0:
            continue
          root_x = _root_x[0][0]
          break

  for i in xrange(len(Hi)):
    if not isinstance(Hi[i], Integer):
      pol_sol = Hi[i].subs(x=root_x, y=yk)
      if not isinstance(pol_sol, Integer):
        roots = pol_sol.roots()
        if len(roots) > 0:
          if len(roots) == 1 and roots[0][0] == 0:
            continue
          _root_y = filter(lambda t: 0 < abs(t[0]) < YY, roots)
          if len(_root_y) == 0:
            continue
          root_y = _root_y[0][0]
          break

  if root_x is None or root_y is None:
    print '[-] Can\'t find solution ... (maybe tweak parameter!)'
    return None, None
  return root_x, root_y

def boneh_durfee_bivariate(_pol, modulo, XX, YY, mm, tt):
  PR = PolynomialRing(ZZ, 'x, y', order='degneglex')
  x, y = PR.gens()
  pol = PR(_pol)
  gik = {}
  hjk = {}
  for k in xrange(mm + 1):
    fk = pol^k * modulo^(mm-k)
    for i in xrange(mm - k + 1):
      gik[i, k] = x^i * fk
    for j in xrange(1, tt + 1):
      hjk[j, k] = y^j * fk

  print gik
  deg_x, deg_y = 0, 0

  monomials = []
  for j in xrange(deg_y):
    for i in xrange(deg_x):
      monomials += [x^i * y^j]
  monomials_xshift = []
  monomials_yshift = []
  for m, monomial in enumerate(monomials):
    if any([gg.monomial_coefficient(monomial) != 0 for gg in gik.values()]):
      monomials_xshift += [monomial]
    elif any([hh.monomial_coefficient(monomial) != 0 for hh in hjk.values()]):
      monomials_yshift += [monomial]
  monomials = monomials_xshift + monomials_yshift
  print monomials

  M = Matrix(ZZ, ((mm + 1) * (mm + 2)) / 2 + (mm + 1) * tt)

  row = 0
  for k in xrange(mm + 1):
    for i in xrange(mm - k + 1):
      for m in xrange(len(monomials)):
        M[row, m] = gik[i, k].monomial_coefficient(monomials[m]) * monomials[m](XX, YY)
      row += 1
  for j in xrange(1, tt + 1):
    for k in xrange(mm + 1):
      for m in xrange(len(monomials)):
        M[row, m] = hjk[j, k].monomial_coefficient(monomials[m]) * monomials[m](XX, YY)
      row += 1
  for k in xrange(mm + 1):
      row += 1

  to_rem = []

  for i in xrange(len(gik) + len(hjk)):
    print M[i]
    if M[i, i] >= modulo ^ mm or all([t == 0 for t in M[i]]):
    # if filter(lambda t: t != 0, M[i])[-1] > modulo ^ mm:
      to_rem += [i]
      pass

  print to_rem
  for ii in to_rem[::-1]:
    M = M[:ii].stack(M[ii+1:])

  matrix_overview(M)

  B = M.LLL()

  matrix_overview(B)

  Hi = {}
  for i in xrange(B.nrows()):
    Hi[i] = 0
    if RR(B[i].norm()) >= modulo/len(monomials):
      continue
    for m in xrange(len(monomials)):
      Hi[i] += ZZ(B[i, m] / monomials[m](XX, YY)) * monomials[m]

  PX = PolynomialRing(ZZ, 'xk')
  PY = PolynomialRing(ZZ, 'yk')
  PK = PolynomialRing(Zmod(modulo), 'xs')

  xk = PX.gen()
  yk = PY.gen()

  root_x = root_y = None

  for i in xrange(len(Hi.keys())):
    for j in xrange(len(Hi.keys())):
      if i == j:
        continue
      pol1 = Hi[Hi.keys()[i]]
      pol2 = Hi[Hi.keys()[j]]
      pol_res = pol1.resultant(pol2, y).subs(x=xk)
      if not isinstance(pol_res, Integer):
        roots = pol_res.roots()
        if len(roots) > 0:
          if len(roots) == 1 and roots[0][0] == 0:
            continue
          _root_x = filter(lambda t: 0 < abs(t[0]) <= XX, roots)
          if len(_root_x) == 0:
            continue
          root_x = (_root_x[0][0]) % modulo
          break

  if root_x is None:
    print '[-] Can\'t find solution for `x0` ... (maybe tweak parameter!)'
    return None, None

  print '[+] `x0` candidate: %s' % root_x
  monic = pol.subs(x=root_x, y=PK.gen()).monic()
  if monic.degree() == 1:
    assert pol(root_x, -monic.constant_coefficient()) % modulo == 0
    return (ZZ(root_x), ZZ(-monic.constant_coefficient()))
  monic = pol.subs(x=root_x, y=PK.gen()).monic()
  print monic.small_roots(X=YY)

  for i in xrange(len(Hi)):
    if not isinstance(Hi[i], Integer):
      pol_sol = Hi[i].subs(x=root_x, y=yk)
      if not isinstance(pol_sol, Integer):
        roots = pol_sol.roots()
        if len(roots) > 0:
          if len(roots) == 1 and roots[0][0] == 0:
            continue
          print roots
          _root_y = filter(lambda t: 0 < abs(t[0]) <= YY, roots)
          if len(_root_y) == 0:
            continue
          # root_y = ZZ(Zmod(modulo)(_root_y[0][0]))
          continue
          break

  if root_x is None or root_y is None:
    print '[-] Can\'t find solution ... (maybe tweak parameter!)'
    return None, None
  return root_x, root_y

def solve_SIP(e, n, delta=0.292, beta=0.5, mm=3):
  F = Zmod(e)
  PR = PolynomialRing(F, 'x, y')
  x, y = PR.gens()
  A = ZZ((n + 1) / 2)
  pol = x * (A + y) - 1
  XX = floor(n ^ delta)
  YY = floor(n ^ beta)
  tt = floor((1 - 2 * delta) * mm)
  x0, y0 = boneh_durfee_bivariate(pol, e, XX, YY, mm, tt)
  if x0 is None or y0 is None:
    return None
  assert pol(x0, y0) == 0
  d0 = (x0 * (A + y0) - 1) / e
  print x0, y0, d0
  assert (Zmod(n)(2)^e)^x0 == 2
  return d0

def main():
  # ASIS CTF 2015 Finals: bodu
  n = 2562256018798982275495595589518163432372017502243601864658538274705537914483947807120783733766118553254101235396521540936164219440561532997119915510314638089613615679231310858594698461124636943528101265406967445593951653796041336078776455339658353436309933716631455967769429086442266084993673779546522240901
  e = 2385330119331689083455211591182934261439999376616463648565178544704114285540523381214630503109888606012730471130911882799269407391377516911847608047728411508873523338260985637241587680601172666919944195740711767256695758337633401530723721692604012809476068197687643054238649174648923555374972384090471828019

  # from wikipedia
  n = 90581
  e = 17993
  print solve_SIP(e, n, delta=0.292, mm=3)


if __name__ == '__main__':
  main()
