from sage.all import *

def factor_shirase_d_3(n):
  '''
  Implementation of yet another Elliptic-Curve Factorization Method proposed by [1] with D = 3.

  References:
  * [1] Masaaki Shirase, 2017, "Condition on composite numbers easily factored with elliptic curve method"
  '''
  F = Zmod(n)
  while True:
    x0 = ZZ.random_element(0, n)
    y0 = ZZ.random_element(0, n)
    B = y0^2 - x0^3
    E = EllipticCurve(F, [0, B])
    # j invariant of `E` is 0
    # because, Hilbert's Class Polynomial for discriminant 3 has root 0.
    P = E(x0, y0)
    try:
      NP = n * P
    except Exception, e:
      # Can't calculate Inverse of some number `k`.
      # iff 1 < g < n then g is non-trivial factor of n where g := gcd(n, k).
      p = gcd(ZZ(e.args[0].split(' ')[2]), n)
      if is_prime(p):
        return p
    x, y, z = map(ZZ, tuple(NP))
    dn2 = gcd(x, y)
    if not dn2.is_square():
      continue
    dn = dn2.nth_root(2)
    if not (x % dn^2 == 0 and y % dn^3 == 0):
      continue
    g = gcd(n, dn)
    if g > 1 and n % g == 0:
        return g


if __name__ == '__main__':
  # HITB AMS 2016 Teaser: Crypto 1000 Special Prime Rib
  n = 0x80fab241e21aacbb5a0c0c58ce7d8a3f844f3f76c2b1006278d79cdd333550ab5f5f86425fdbf06063481d7d7922f1c17083532285b1d8faee843d8a02e74f277a47084bc5585f0d16a2ab7f2e2c074a274c9b890b05a4ed05739f9baeaa501c265d68c04c146a5daed6ef5e0a45aa7c9ae1e7c3741c39f7f00936d1d627bc5b
  p = factor_shirase_d_3(n)
  q = n / p
  print p, q
  assert n == p*q
