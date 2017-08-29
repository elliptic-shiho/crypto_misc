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


def coron_bivariate_integer_small_root(poly, XX, YY, kk):
  '''
  Solve Bivariate Polynomial Small Root

  Implementation of [1].

  References:
    * [1] Jean-Se'bastien Coron. 2004. "Finding Small Roots of Bivariate Integer Polynomial Equations Revisited"

  Author:
    Hayato Ashida
'''

  p00 = poly.constant_coefficient()
  assert gcd(p00, XX * YY) == 1, 'p00 and XY has common divisor'
  x, y = poly.parent().gens()

  monomials = list(poly.monomials())
  monomials.sort()

  WW = max(map(lambda t: abs(poly.monomial_coefficient(t)) * t(XX, YY), monomials))
  uu = WW + 1 + ZZ((1 - WW) % abs(p00))
  delta = max(poly.degree(x), poly.degree(y))
  omega = (delta + kk + 1)^2
  nn = uu * (XX * YY)^kk

  print '[+] Bound Check...',
  sys.stdout.flush()
  if RR(XX*YY) < RR(WW^((2/3) * delta)):
    print 'OK'
  else:
    print 'Failed (maybe not found solution...)'

  if p00 != 0:
    F = Zmod(nn)
    PF = PolynomialRing(F, 'xn, yn')
    q = poly.parent()(PF(poly) * F(p00)^-1)

  # Construct Polynomial for lattice construction (cf. [1] p.7)
  qij = {}
  for i in xrange(0, kk + 1):
    for j in xrange(0, kk + 1):
      qij[i, j] = x^i * y^j * XX^(kk-i) * YY^(kk-j) * q
  index_range = sorted(list(set(itertools.product(range(0, delta + kk + 1), repeat=2)) - set(itertools.product(range(0, kk + 1), repeat=2))))
  for i, j in index_range:
    qij[i, j] = x^i * y^j * nn

  monomials = set()
  for k in qij.keys():
    monomials |= set(qij[k].monomials())
  monomials = list(monomials)
  monomials.sort()

  M = Matrix(ZZ, omega, omega)
  assert len(monomials) == omega

  # Construct Lattice
  col = 0
  for i in xrange(kk + 1):
    for j in xrange(kk + 1):
      q_cur = qij[i, j]
      for ii, m in enumerate(monomials):
        M[col, ii] = q_cur.monomial_coefficient(m) * m(XX, YY)
      col += 1

  for i, j in index_range:
    q_cur = qij[i, j]
    for ii, m in enumerate(monomials):
      M[col, ii] = q_cur.monomial_coefficient(m) * m(XX, YY)
    col += 1

  matrix_overview(M)

  print '\n===\n'

  # LLL
  B = M.LLL()

  matrix_overview(B)
  print

  # Solve equatation for each variable
  PK = PolynomialRing(ZZ, 'xk, yk')
  xk, yk = PK.gens()

  PX = PolynomialRing(ZZ, 'xs')
  xs = PX.gen()
  PY = PolynomialRing(ZZ, 'ys')
  ys = PY.gen()

  monomials = map(lambda t: PK(t), monomials)
  pkf = PK(poly)
  x_root = y_root = None

  # Re-construct polynomial from LLL-reduced matrix `B`
  H = [(i, 0) for i in xrange(omega)]
  H = dict(H)
  for i in xrange(omega):
    for j in xrange(omega):
      H[i] += PK((monomials[j] * B[i, j]) / monomials[j](XX, YY))

  # Solve for `x`
  # My Heuristics: finding resultant from all polynomials
  for i in xrange(omega):
    pol = H[i].resultant(pkf, yk).subs(xk=xs)
    if not isinstance(pol, Integer):
      roots = pol.roots()
      roots = filter(lambda t: 0 < t[0] < XX, roots)
      if len(roots) != 0:
        print '[+] Found Solution for x'
        x_root = roots[0][0]
        break
  else:
    print '[+] solution not found...'
    return None, None

  # Solve for `y`
  roots_y = PY(pkf.subs(xk=x_root).subs(yk=ys)).roots()
  roots_y = filter(lambda t: 0 < t[0] < YY, roots_y)
  if len(roots_y) != 0:
    print '[+] Found Solution for y'
    y_root = roots_y[0][0]
  return x_root, y_root


def main():
  # Sample Implementation: High-bit known attack for RSA (cf. [1] Appendix C)
  p = next_prime(ZZ.random_element(2^255, 2^256))
  q = next_prime(ZZ.random_element(2^255, 2^256))
  n = p*q
  x0 = p % 2^64
  y0 = q % 2^64
  p0 = p - x0
  q0 = q - y0

  XX = next_prime(2^64)
  YY = next_prime(2^64)
  kk = 1

  PR = PolynomialRing(ZZ, 'x, y')
  x, y = PR.gens()

  f = (p0 * q0 - n) + q0 * x + p0 * y + x*y
  assert f(x0, y0) == 0
  roots = coron_bivariate_integer_small_root(f, XX, YY, kk)
  assert f(roots[0], roots[1]) == 0


if __name__ == '__main__':
  main()


'''
Fri Apr 28 01:56:19 JST 2017 ~/prog/lab/cryptography/crypto-misc/small_root 100%
> time sage coron.sage
[+] Bound Check... OK
00 X X X   X         
01   X   X X   X     
02     X   X X   X   
03         X   X X X 
04       X           
05             X     
06           X       
07               X   
08                 X 

===

00 X X X X X X X X X 
01 X X X X X X X X X 
02 X X X X X X X X X 
03 X X X X X X X X X 
04 X X X X X X X X X 
05 X X X X X X X X X 
06 X X X X X X X X X 
07 X X X X X X X X X 
08 X X X X X X X X X 

[+] Found Solution for x
[+] Found Solution for y

real    0m4.509s
user    0m4.308s
sys     0m0.136s
'''
