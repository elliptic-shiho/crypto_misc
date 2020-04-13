from sage.all import *

# Qn_tau is implementation of Quotient ring Z/nZ[x] / (x^2 - tau)
load('Qn_tau.sage')
# EC is simple implementation of Elliptic Curve
# (because, sage's elliptic curve implementation is not compatible with Qn_tau)
load('EC.sage')

def factor_shirase_linear(n, D):
  '''
  Implementation of yet another Elliptic-Curve Factorization Method proposed by [1]
    with D = 11, 19, 43, 67, 163, if and only if H_D(j) is linear
      where H_D(j) is Hilbert's Class Polynomial with discriminant is D.

  References:
  * [1] Masaaki Shirase, 2017, "Condition on composite numbers easily factored with elliptic curve method"
  '''
  F = Zmod(n)
  PF = PolynomialRing(F, 'X')
  x = PF.gen()
  H_D = hilbert_class_polynomial(-D)
  assert H_D.degree() == 1, 'Class polynomial must have a degree 1'
  j0 = H_D.roots()[0][0]
  j0_inv_1728 = ZZ(inverse_mod(1728 - j0, n))
  print '[+] j0 = %d' % j0
  while True:
    R = ZZ.random_element(1, n)
    ADR = (3 * j0 * R^2 * j0_inv_1728) % n
    BDR = (2 * j0 * R^3 * j0_inv_1728) % n
    x0 = ZZ.random_element(0, n)
    tau = (x0^3 + ADR * x0 + BDR) % n
    FQ = Qn_tau(n, tau)
    E = EC(FQ, ADR, BDR)
    P = ECPoint(FQ(x0, 0), FQ(0, 1), 1)
    try:
      nP = E.mul(n, P)
    except ZeroDivisionError, e:
      t = e.args[0].split()
      x0, x1 = ZZ(t[0]), ZZ(t[2].replace('X', ''))
      g = gcd(x0^2 - x1^2 * tau, n)
      if 0 < g < n:
        return g

if __name__ == '__main__':
  # HITB AMS 2016 Teaser: Crypto 1000 Special Prime Rib
  n = 0x80fab241e21aacbb5a0c0c58ce7d8a3f844f3f76c2b1006278d79cdd333550ab5f5f86425fdbf06063481d7d7922f1c17083532285b1d8faee843d8a02e74f277a47084bc5585f0d16a2ab7f2e2c074a274c9b890b05a4ed05739f9baeaa501c265d68c04c146a5daed6ef5e0a45aa7c9ae1e7c3741c39f7f00936d1d627bc5b
  n = 122885643723000432249644760468389188624901340133395238224842000337736669981443704071141388121472030756208165869797637031888222362809900683103944835376737386357640106122079313347291918461873652446934855347009443865079947354999865318372362290819783827546603774471438105605546446537286320695311365087282613747067
  p = factor_shirase_linear(n, 11)
  assert n % p == 0
  assert 0 < p < n
  q = n / p
  assert p * q == n
