from sage.all import *

def matrix_overview(BB):
  # display matrix picture with space and X
  # references: https://github.com/mimoo/RSA-and-LLL-attacks/blob/master/boneh_durfee.sage
  for ii in range(BB.dimensions()[0]):
    a = ('%02d ' % ii)
    for jj in range(BB.dimensions()[1]):
      a += ' ' if BB[ii,jj] == 0 else 'X'
      if BB.dimensions()[0] < 60:
        a += ' '
    print a

def howgrave_graham_univariate(_pol, modulo, XX, hh):
  '''
  Solve Modular Univariate Polynomial: Howgrave-Graham's Method

  References:

  * [1] Nicholas Howgrave-Graham. 1997. "Finding small roots of univariate modular equations revisited"
  '''
  assert hh >= 2
  assert _pol.is_monic(), 'polynomial must be monic'
  PK = PolynomialRing(ZZ, 'x')
  x = PK.gen()
  pol = PK(_pol)

  kk = pol.degree()
  M = Matrix(ZZ, hh * kk)

  quv = {}
  for i in xrange(hh * kk):
    v = floor(RR(i/kk))
    u = i - kk * v
    quv[u, v] = modulo ^ (hh - 1 - v) * x^u * pol^v

  for i in xrange(hh * kk):
    for j in xrange(hh * kk):
      v = floor(RR(i/kk))
      u = i - kk * v
      M[i, j] = quv[u, v].monomial_coefficient(x^j) * XX ^ j

  print '[+] Matrix M:'
  matrix_overview(M)
  print

  B = M.LLL()

  print '[+] Matrix B:'
  matrix_overview(B)
  print 

  pol = PK(reduce(lambda u, t: u + (x^t[0] * t[1]) / (XX^t[0]), enumerate(B[0]), 0))
  roots = pol.roots()

  assert len(roots) > 0, 'Can\'t Find Solution!!'
  if len(roots) == 1 and roots[0][0] == 0:
    print 'couldn\'t find non-trivial solution...'
    return [0]
  return map(lambda t: t[0], roots)

def main():
  modulo = 35
  F = Zmod(modulo)
  PR = PolynomialRing(F, 'x')
  x = PR.gen()
  pol = x^2 + 14*x + 19

  assert pol(3) == 0

  XX = 4
  roots = howgrave_graham_univariate(pol, modulo, XX, 3)
  assert len(roots) > 0
  x0 = roots[0]
  assert pol(x0) == 0

if __name__ == '__main__':
  main()
