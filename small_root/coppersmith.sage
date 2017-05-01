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

def coppersmith_univariate(_pol, modulo, XX, hh):
  '''
  Implementaion of Coppersmith's Solving Univariate Modular Equation with Small Root [1].

  References:

    [1] Don Coppersmith. 1996. "Finding a Small Root of a Univariate Modular Equation"
  '''
  _pol = _pol.monic()
  kk = _pol.degree()
  epsilon = 0.1
  assert hh >= floor(max(RR(7/kk), RR((kk + epsilon * kk - 1) / (epsilon * kk^2))))
  PR = PolynomialRing(ZZ, 'x')
  x = PR.gen()
  pol = _pol.change_ring(ZZ)
  qij = {}
  for i in xrange(0, kk):
    for j in xrange(1, hh):
      qij[i, j] = x^i * pol^j

  A = Matrix(QQ, hh * kk, hh * kk - kk)

  for g in xrange(hh * kk):
    for i in xrange(kk):
      for j in xrange(1, hh):
        '''
        In the paper, `gamma(i, j)` is defined by hk+i+(j-1)k.
        Actually, index of a matrix needs range 0 <= `gamma(i, j)` <= hk".
        However, `gamma(i, j)` greater than `hk`. so, we set
        `gamma'(i, j) := gamma(i, j) - hk` to match a index of matrix and use it.
        (cf. [1] p.157)
        '''
        A[g, i + (j - 1) * kk] = qij[i, j].monomial_coefficient(x^g)

  D_ = Matrix(QQ, hh * kk - kk)
  for i in xrange(hh * kk - kk):
      D_[i, i] = modulo^(floor((i+1)/hh) + 1)
  print D_

  D = Matrix(QQ, hh * kk)
  delta = (RR(sqrt(hh * kk)) * 10000).integer_part() / 10000
  delta = 1
  for g in xrange(hh * kk):
    D[g, g] = delta * (1 / (XX^g))

  print D
  # Zero Matrix hk x hk
  O = Matrix(QQ, (hh - 1) * kk, hh * kk)

  M = (D.stack(O)).augment(A.stack(D_))

  matrix_overview(M)

  assert M.ncols() == (2*hh*kk - kk)
  assert M.nrows() == (2*hh*kk - kk)

  HM = copy(M)
  # Swap some rows: prepare for identity submatrix - lower right
  for i in xrange(0, hh * kk):
    for g in xrange(2*hh*kk - kk - 1, hh*kk-1, -1):
      HM.swap_rows(i, g)

  # Do elementary row operations (like Gaussian Elimination)
  # for transform to identity matrix at lower-right.
  for i in xrange(hh*kk, 2*hh*kk - kk):
    for g in xrange(hh*kk, 2*hh*kk - kk):
      pivot = HM[g, i]
      if pivot != 0:
        if g == i:
          pivot = pivot - 1
        HM.add_multiple_of_row(g, i, -pivot)

  # Do elementary row operations (like Gaussian Elimination)
  # for transform to zero matrix at lower-left.
  for i in xrange(hh*kk, 2*hh*kk - kk):
    for g in xrange(0, hh * kk):
      pivot = HM[g, i]
      if pivot != 0:
        HM.add_multiple_of_row(g, i, -pivot)

  print '[+] Finish transform.'
  matrix_overview(HM)

  M_hat = HM[0:hh * kk, 0:hh * kk]

  # transform to triangle matrix
  M_hat = M_hat[-kk:].stack(M_hat[:-kk])

  print '[+] \hat{M}'
  matrix_overview(M_hat)

  # LLL
  B = M_hat.LLL()

  print B

if __name__ == '__main__':
  N = 35
  F = Zmod(N)

  PR = PolynomialRing(F, 'x')
  x = PR.gen()

  pol = x^3 + 14*x + 19
  print pol
  coppersmith_univariate(pol, N, 2, 3)

