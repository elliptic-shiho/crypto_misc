from sage.all import *

def factor_shirase_linear(n, D):
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
    Q_tau = PF.quo(x^2 - tau)
    E = EllipticCurve(F, [ADR, BDR])
    EQ = E.base_extend(Q_tau)
    print EQ
    X = Q_tau.gen()
    P = EQ(x0, X)
    print X^2 - (x0^3 + ADR * x0 + BDR)
    eq = pari(EQ)
    print eq


if __name__ == '__main__':
  # HITB AMS 2016 Teaser: Crypto 1000 Special Prime Rib
  n = 0x80fab241e21aacbb5a0c0c58ce7d8a3f844f3f76c2b1006278d79cdd333550ab5f5f86425fdbf06063481d7d7922f1c17083532285b1d8faee843d8a02e74f277a47084bc5585f0d16a2ab7f2e2c074a274c9b890b05a4ed05739f9baeaa501c265d68c04c146a5daed6ef5e0a45aa7c9ae1e7c3741c39f7f00936d1d627bc5b
  factor_shirase_linear(n, 11)
