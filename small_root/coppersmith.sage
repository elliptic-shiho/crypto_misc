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

def coppersmith_univariate(_pol, modulo, XX=None, hh=None, epsilon=None):
  '''
  Implementaion of Coppersmith's Solving Univariate Modular Equation with Small Root [1].

  References:

    [1] Don Coppersmith. 1996. "Finding a Small Root of a Univariate Modular Equation"
    [2] Nicholas Howgrave-Graham. 1997. "Finding Small Roots of Univariate Modular Equations Revisited"
  '''
  assert _pol.is_monic(), 'p(x) must be monic'
  kk = _pol.degree()
  if epsilon is None:
    epsilon = RR(1/log(modulo, 2))
  if hh is None:
    hh = floor(max(RR(7/kk), RR((kk + epsilon * kk - 1) / (epsilon * kk^2)))) + 1
  if XX is None:
    XX = floor(RR(0.5 * modulo ^ (1/k - epsilon)))
  assert hh >= floor(max(RR(7/kk), RR((kk + epsilon * kk - 1) / (epsilon * kk^2))))
  PR = PolynomialRing(ZZ, 'x')
  x = PR.gen()
  pol = PR(_pol)
  qij = {}
  for i in xrange(0, kk):
    for j in xrange(1, hh):
      qij[i, j] = x^i * pol^j

  # Matrix Construction: Matrix of coefficient-vectors
  A = Matrix(QQ, hh * kk, hh * kk - kk)

  for g in xrange(hh * kk):
    for i in xrange(kk):
      for j in xrange(1, hh):
        '''
        In the paper[1], `gamma(i, j)` is defined by hk+i+(j-1)k.
        Actually, index of a matrix needs range 0 <= `gamma(i, j)` <= hk".
        However, `gamma(i, j)` greater than `hk`. so, we set
        `gamma'(i, j) := gamma(i, j) - hk` to match a index of matrix and use it.
        (cf. [1] p.157)
        '''
        A[g, i + (j - 1) * kk] = qij[i, j].monomial_coefficient(x^g)

  # Matrix Construction: Diagonal Matrix. D' = (d'_{ij}) (0 <= i,j < (h-1)k) where d'_{ij} = N^{floor((i + 1) / h) + 1} if i = j, otherwise d'_{ij} = 0. 
  D_ = Matrix(QQ, hh * kk - kk)
  for i in xrange(hh * kk - kk):
      D_[i, i] = modulo^(floor((i+1)/hh) + 1)

  # Matrix Construction: Diagonal Matrix. D = (d_{ij}) (0 <= i, j < hk) where d_{ij} = delta * (X^-i) if i = j, otherwise d_{ij} = 0

  D = Matrix(QQ, hh * kk)
  delta = QQ((RR(sqrt(hh * kk)) * 1000000).integer_part() / 1000000)
  # delta = 1

  for g in xrange(hh * kk):
    D[g, g] = delta * (1 / (XX^g))

  # Zero Matrix (h-1)k x hk
  O = Matrix(QQ, (hh - 1) * kk, hh * kk)

  '''
  Combined Matrix (2h-1)k x (2h-1)k

      |           |
      |  D  |  A  |
  M = | ----+---- |
      |  O  |  D' |
      |           |
  '''

  M = (D.stack(O)).augment(A.stack(D_))

  print '[+] Matrix M'
  matrix_overview(M)

  assert M.ncols() == (2*hh*kk - kk)
  assert M.nrows() == (2*hh*kk - kk)

  # Transform Matrix to extract sub-matrix \hat{M} (|det(M)| = |det(\hat{M})|)
  HM = copy(M)

  # Swap some rows: prepare for identity submatrix - lower right
  for i in xrange(0, hh * kk):
    for g in xrange(2*hh*kk - kk - 1, hh*kk-1, -1):
      HM.swap_rows(i, g)

  # Swapping rows `0` to `hk-k` and `hk-k` to `hk`.
  HM = HM[hh*kk-kk:hh*kk].stack(HM[0:hh*kk-kk]).stack(HM[hh*kk:])

  '''
  Do elementary row operations (like Gaussian Elimination)

  for transform to identity matrix at lower-right.
  Note: p(x) is monic, we know all q_{ij} is monic. so, all diagonal element of upper-right submatrix `A` is also 1.

        |           |
        |  *  |  *  |
  M_1 = | ----+---- |
        |  A' |  I  |
        |           |
  '''

  for i in xrange(hh*kk, 2*hh*kk - kk):
    for g in xrange(hh*kk, 2*hh*kk - kk):
      pivot = HM[g, i]
      if pivot != 0:
        if g == i:
          pivot = pivot - 1
        HM.add_multiple_of_row(g, i, -pivot)

  '''
  Do elementary row operations for transform to zero matrix at upper-right.

        |           |
        |  M^ |  O  |
  M_2 = | ----+---- |
        |  A' |  I  |
        |           |
  '''

  for i in xrange(hh*kk, 2*hh*kk - kk):
    for g in xrange(0, hh * kk):
      pivot = HM[g, i]
      if pivot != 0:
        HM.add_multiple_of_row(g, i, -pivot)

  print '[+] Transformed Matrix:'
  matrix_overview(HM)

  # \hat{M} (i.e. M^)
  M_hat = HM[0:hh * kk, 0:hh * kk]

  # transform to lower-triangle matrix (cf. [2] p.137)
  M_hat, denom = (M_hat)._clear_denom()
  for i in xrange(hh*kk):
    M_hat[i] = M_hat[i][::-1]
  M_hat = M_hat[::-1]

  print '[+] \hat{M}:'
  matrix_overview(M_hat)

  # LLL
  B = M_hat.LLL()

  print '[+] B (LLL-reduced matrix):'
  matrix_overview(B)

  # Create orthogonal matrix for checking bounds (cf. [1] p.158, [2] p.132)
  Bstar_vec = [B[0]]
  for i in xrange(1, B.nrows()):
    s = B[i]
    for j in xrange(0, i):
      s -= (Bstar_vec[j].dot_product(B[i]) / Bstar_vec[j].dot_product(Bstar_vec[j])) * Bstar_vec[j]
    Bstar_vec += [s]

  # Bound Check (This condition is always True)
  assert RR(Bstar_vec[-1].norm()) > RR(M_hat.det() ^ (1/(hh*kk)) * 2 ^ (-((hh*kk)-1)/4)) > 1

  H_1 = M * HM^-1
  # H_2 = M_hat * Matrix(Bstar_vec)^-1
  H_2 = M_hat * B^-1

  PQ = PolynomialRing(QQ, 'xq')
  xq = PQ.gen()
  pol_qq = PQ(pol)

  # Constuct Polynomial Vector: (1, x0, ..., x0^{hk-1}, -y0, -y0x0, ..., -y0x0^{k-1}, ..., -y0^{h-1}x0^{k-1})
  # cf. [1] pp.157 - 159, [2] pp.135 - 137
  v = []
  v = Matrix(PQ, 1, 2*hh*kk-kk)
  for i in xrange(hh * kk):
    v[0, i] = xq^i

  for i in xrange(1, hh):
    for j in xrange(kk):
      v[0, hh*kk + (i - 1)*kk + j] = -xq^j * (pol_qq^i / modulo^i)

  # Zero check (This condition is always True)
  assert all([t == 0 for t in (v * H_1)[0][hh * kk:]])

  pol_sol = (v * H_1)[0][:hh*kk] * (H_2.T)[hh * kk - 1][::-1]

  # Finding roots for `x`
  roots = []
  for x_candidate in map(lambda t: t[0], (pol_sol).roots()):
    if x_candidate.denom() != 1:
      continue
    x_candidate = ZZ(x_candidate)
    if _pol(x_candidate) == 0:
      roots += [x_candidate]
  return roots

def main():
  '''
  [2] 4.1 Examples: Coppersmith's Method

  The equation `x^2 + 14x + 19 = 0 mod 35` has small root `x0 = 3`. (this is small root!)
  using coppersmith's method, we can solve this equation.
  '''
  N = 35
  F = Zmod(N)

  PR = PolynomialRing(F, 'x')
  x = PR.gen()

  pol = x^2 + 14*x + 19
  print pol

  assert coppersmith_univariate(pol, N, 4, 3) == [3]

if __name__ == '__main__':
  main()

'''
Wed May  3 04:18:09 JST 2017 ~/prog/lab/cryptography/crypto-misc/small_root 100%
> time sage coppersmith.sage
x^2 + 14*x + 19
[+] Matrix M
00 X           X   X   
01   X         X X X X 
02     X       X X X X 
03       X       X X X 
04         X       X X 
05           X       X 
06             X       
07               X     
08                 X   
09                   X 

[+] Transformed Matrix:
00 X   X X X X         
01   X X X X X         
02     X X X X         
03       X X X         
04         X X         
05           X         
06     X X X X X       
07       X X X   X     
08         X X     X   
09           X       X 

[+] \hat{M}:
00 X           
01 X X         
02 X X X       
03 X X X X     
04 X X X X X   
05 X X X X   X 

[+] B (LLL-reduced matrix):
00 X X X X X   
01 X X X X X   
02 X X X       
03 X           
04 X X X X X X 
05 X X X X X X 


real    0m4.627s
user    0m4.360s
sys     0m0.208s
'''

