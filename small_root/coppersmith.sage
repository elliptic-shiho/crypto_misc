from sage.all import *
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
  kk = _pol.degree()
  epsilon = 0.1
  assert hh >= floor(max(RR(7/kk), RR((kk + epsilon * kk - 1) / (epsilon * kk^2))))
  PR = PolynomialRing(ZZ, 'x')
  x = PR.gen()
  pol = _pol.change_ring(ZZ)
  qij = {}
  for i in xrange(0, kk + 1):
    for j in xrange(1, hh):
      qij[i, j] = x^i * pol^j

  D = Matrix(QQ, hh*kk)

  for i in xrange(hh*kk):
    D[i, i] = 1/XX^i

  O = Matrix(QQ, (hh-1)*kk, hh*kk)

  D_ = Matrix(QQ, (hh-1) * kk)

  for i in xrange((hh-1)*kk):
    D_[i, i] = modulo^(floor((kk + i) / kk))

  A = Matrix(QQ, hh * kk, (hh - 1) * kk)
  pj = {}
  for j in xrange(1, (hh - 1) * kk + 1):
    v = floor((kk + j - 1) / kk)
    u = (j - 1) - kk * (v - 1)
    pj[j - 1] = (x^u * pol^v)

  for i in xrange(hh * kk):
    for j in xrange((hh - 1) * kk):
      A[i, j] = pj[j].monomial_coefficient(x^(i))

  M = (D.stack(O)).augment(A.stack(D_))

  matrix_overview(M)

if __name__ == '__main__':
  N = 35
  F = Zmod(N)

  PR = PolynomialRing(F, 'x')
  x = PR.gen()

  pol = x^2 + 14*x + 19
  print pol
  coppersmith_univariate(pol, N, 2, 3)

