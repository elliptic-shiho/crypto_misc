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

def boneh_durfee_bivariate(_pol, modulo, mm, tt):
  PR = PolynomialRing(ZZ, 'x, y', order='degneglex')
  x, y = PR.gens()
  pol = PR(_pol)
  gik = {}
  hjk = {}
  for k in xrange(mm + 1):
    fk = pol^k * modulo^(mm-k)
    for i in xrange(mm - k):
      gik[i, k] = x^i * fk
    for j in xrange(tt + 1):
      hjk[j, k] = y^j * fk

  monomials = set()
  for k in gik.keys():
    monomials |= set(gik[k].monomials())
  for k in hjk.keys():
    monomials |= set(hjk[k].monomials())
  monomials = list(monomials)
  monomials.sort()
  print monomials
  print len(monomials)
  M = Matrix(ZZ, len(monomials))


def solve_SIP(e, n):
  F = Zmod(e)
  PR = PolynomialRing(F, 'x, y')
  x, y = PR.gens()
  A = ZZ((n + 1) / 2)
  pol = (x * (A + y) - 1)
  return boneh_durfee_bivariate(pol, e, 3, 1)

def main():
  # ASIS CTF 2015 Finals: bodu
  n = 2562256018798982275495595589518163432372017502243601864658538274705537914483947807120783733766118553254101235396521540936164219440561532997119915510314638089613615679231310858594698461124636943528101265406967445593951653796041336078776455339658353436309933716631455967769429086442266084993673779546522240901
  e = 2385330119331689083455211591182934261439999376616463648565178544704114285540523381214630503109888606012730471130911882799269407391377516911847608047728411508873523338260985637241587680601172666919944195740711767256695758337633401530723721692604012809476068197687643054238649174648923555374972384090471828019
  print solve_SIP(e, n)


if __name__ == '__main__':
  main()
