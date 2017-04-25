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

def solve_modular_univariate_polynomial(pol, modulo, XX):
  dd = pol.degree()
  PK = PolynomialRing(ZZ, 'xk')
  f_ = (pol * pol[0]^-1)
  f_ = f_.change_ring(PK)

  coeffs = list(f_)
  M = Matrix(ZZ, dd + 1)

  monomials = sorted(f_.monomials())

  for i in xrange(dd + 1):
    M[0, i] = coeffs[i] * monomials[i](XX)

  for i in xrange(1, dd + 1):
    M[i, i] = modulo * monomials[i](XX)

  print '[+] Matrix M:'
  matrix_overview(M)
  print

  B = M.LLL()

  print '[+] Matrix B:'
  matrix_overview(B)
  print 

  pol = PK(reduce(lambda u, t: u + (t[0] * t[1]) / t[0](XX), zip(monomials, B[0]), 0))
  roots = pol.roots()

  assert len(roots) > 0, 'Can\'t Find Solution!!'
  if len(roots) == 1 and roots[0][0] == 0:
    print 'couldn\'t find non-trivial solution...'
    return [0]
  return map(lambda t: t[0], roots)

def main():
  p = next_prime(2^512)
  q = next_prime(p)
  n = p*q

  beta = 0.2
  dd = 2

  root = ZZ.random_element(1, floor(n^beta))
  print '[+] root = %d' % root

  F = Zmod(n)
  PR = PolynomialRing(F, 'x')
  x = PR.gen()
  pol = expand(PR.random_element(degree=(-1, dd-1)) * (x-root))
  print '[+] Polynomial = %s' % pol
  assert pol(root) == 0

  XX = floor(n^beta)
  roots = solve_modular_univariate_polynomial(pol, n, XX)
  assert len(roots) > 0
  assert root in roots
  print '[+] check OK'

if __name__ == '__main__':
  main()
