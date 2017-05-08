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

def is_geometrically_progressive(_M, a, b, C, D, c0, c1, c2, c3, c4, beta):
  assert _M.ncols() == (a + 1) * b
  assert _M.nrows() == (a + 1) * b
  M = lambda i,j,k,l: _M[b * i + j, b * k + l]
  if not (beta * c1 + c3 >= 0 and beta * (c2 + c4)):
    return False
  for i in xrange(a + 1):
    for j in xrange(1, b + 1):
      for k in xrange(a + 1):
        for l in xrange(1, b + 1):
          if i == k and j == l:
            if not M(k, l, k, l) == D ^ (c0 + c1 * k + c2 * l + c3 * k + c4 * l):
              return False
          elif i > k and j > l:
            if not M(i, j, k, l) == 0:
              return False
          else:
            if not abs(M(i, j, k, l)) <= C * D ^ (c0 + c1 * i + c2 * j + c3 * k + c4 * l):
              return False

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
  PR = PolynomialRing(ZZ, 'x, y', order='deglex')
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

  monomials_xshift = set()
  for k in gik.keys():
    monomials_xshift |= set(gik[k].monomials())
  monomials_yshift = set()
  for k in hjk.keys():
    monomials_yshift |= set(hjk[k].monomials())
  monomials_yshift -= monomials_xshift
  monomials_xshift = list(monomials_xshift)
  monomials_xshift.sort()
  monomials_yshift = list(monomials_yshift)
  monomials_yshift.sort()

  monomials = monomials_xshift + monomials_yshift
  print monomials

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

  to_rem = []

  for i in xrange(len(hjk)):
    ii = len(gik) + i
    # if M[ii, ii] >= modulo ^ mm:
    if filter(lambda t: t != 0, M[ii])[-1] > modulo ^ mm:
      to_rem += [ii]

  print to_rem
  #for ii in to_rem[::-1]:
    #M = M[:ii].stack(M[ii+1:])

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
          _root_x = filter(lambda t: 0 < abs(t[0]) <= XX, roots)
          if len(_root_x) == 0:
            continue
          root_x = _root_x[0][0]
          break

  if root_x is None:
    print '[-] Can\'t find solution for `x0` ... (maybe tweak parameter!)'
    return None, None

  print '[+] `x0` candidate: %s' % root_x

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
          root_y = _root_y[0][0]
          break

  if root_x is None or root_y is None:
    print '[-] Can\'t find solution ... (maybe tweak parameter!)'
    return None, None
  return root_x, root_y

def solve_SIP(e, n, delta=0.3, beta=0.5, mm=5):
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
  assert (Zmod(n)(2)^e)^d0 == 2
  return d0

def main():
  # ASIS CTF 2015 Finals: bodu
  n = 2562256018798982275495595589518163432372017502243601864658538274705537914483947807120783733766118553254101235396521540936164219440561532997119915510314638089613615679231310858594698461124636943528101265406967445593951653796041336078776455339658353436309933716631455967769429086442266084993673779546522240901
  e = 2385330119331689083455211591182934261439999376616463648565178544704114285540523381214630503109888606012730471130911882799269407391377516911847608047728411508873523338260985637241587680601172666919944195740711767256695758337633401530723721692604012809476068197687643054238649174648923555374972384090471828019

  # from wikipedia
  #n = 90581
  #e = 17993
  print solve_SIP(e, n)


if __name__ == '__main__':
  main()
