'''
An implementation of [1]. 

[1]: Matus Nemec, Marek Sys, Patr Svenda, Dusan Klinec, and Vashek Matyas. 2017. "The Return of Coppersmith's Attack: Practical Factorization of Widely Used RSA Moduli"
[2]: Alexander May. 2003. "New RSA Vulnerabilities UsingLattice Reduction Methods"
'''
import functools
import binascii
import math

# display matrix picture with 0 and X
# references: https://github.com/mimoo/RSA-and-LLL-attacks/blob/master/boneh_durfee.sage
def matrix_overview(BB):
  for ii in range(BB.dimensions()[0]):
    a = ('%02d ' % ii)
    for jj in range(BB.dimensions()[1]):
      a += ' ' if BB[ii,jj] == 0 else 'X'
      if BB.dimensions()[0] < 60:
        a += ' '
    print(a)
  print("")

def factor_expand(n):
  res = []
  for p, e in factor(n):
    for i in range(0, e):
      res += [p**(i + 1)]
  return res

def reward_at_cost(old_M, new_M, old_order, new_order):
  '''
  Cost function
  cf. [1] 2.7.2 Greedy heuristic
  '''
  return (log(old_order, 2) - log(new_order, 2)) / (log(old_M, 2) - log(new_M, 2))

def compute_new_M_from_order(M, order_new_M):
  '''
  Compute M' from M, such that ord_{M'}(65537) = order_new_M. 
  cf. [1] Algorithm 2. 
  '''
  M_ = M
  for pi in factor_expand(M):
    order_pi = Mod(65537, pi).multiplicative_order()
    if order_new_M % order_pi != 0:
      M_ = M_ // pi
  return M_

def find_new_M(M, log2N):
  '''
  Find M' from M
  '''
  M_ = M
  order = Mod(65537, M_).multiplicative_order()
  while True:
    L = []
    for fi in factor_expand(order):
      new_order = order // fi
      new_M = compute_new_M_from_order(M_, new_order)
      L += [(reward_at_cost(M_, new_M, order, new_order), (new_order, new_M))]
    _, (new_order, new_M) = list(reversed(sorted(L)))[0]
    if log(new_M, 2) < log2N / 4:
      break
    M_ = new_M
    order = new_order
  return M_

def coppersmith_univariate(pol, NN, XX, mm, tt, beta=1.0):
  '''
  An implementation of Coppersmith's method for univariate polynomial

  Note: I didn't use [1]'s definition. This implementation referenced [2]. 
  '''
  PR.<x> = PolynomialRing(ZZ)
  polZZ = PR(pol)
  delta = polZZ.degree()
  pols = []

  for i in range(mm):
    for j in range(delta):
      pols += [(x * XX)^j * NN^(mm - i) * polZZ(x * XX)^i]

  for i in range(tt):
    pols += [(x * XX)^i * polZZ(x * XX)^mm]

  deg = delta * mm + tt
  M = Matrix(ZZ, deg)
  for i in range(deg):
    for j in range(deg):
      M[i, j] = pols[i].monomial_coefficient(x^j)

  B = M.LLL()

  f = 0
  for i in range(deg):
    f += x^i * B[0, i] / XX^i

  roots = []
  for x0, _ in f.roots():
    if x0.is_integer():
      if gcd(polZZ(ZZ(x0)), NN) >= NN^beta:
        roots += [ZZ(x0)]
  return list(set(roots))

def roca_attack(N, M_, mm, tt):
  '''
  ROCA Attack

  * mm and tt are tweakable parameter

  cf. [1] Algorithm 1.
  '''
  c_ = discrete_log(N, Mod(65537, M_))
  order_ = Mod(65537, M_).multiplicative_order()
  a_ = c_ // 2
  upper_bounds = (c_ + order_) // 2
  ZmodN = Zmod(N)
  PR.<x> = PolynomialRing(ZmodN, implementation='NTL')
  print("[+] c' = {}".format(c_))
  print("[+] Iteration range: [{}, {}]".format(a_, upper_bounds))
  while a_ <= upper_bounds:
    const = (Mod(65537, M_)^a_).lift()
    f = x + Mod(M_, N)^-1 * const
    beta = 0.5
    XX = floor(2 * N^beta / M_)
    # small_roots is useless (It can't solve this polynomial). 
    # res = f.small_roots(beta=beta, X=XX) 
    # Pari/GP's zncoppersmith function can solve
    res = coppersmith_univariate(f, N, XX, mm, tt, beta)
    for k_ in res:
      p = k_ * M_ + const
      if N % p == 0:
        return (p, N // p)
    a_ += 1

def main():
  '''
  Test
  '''
  M = 0x924cba6ae99dfa084537facc54948df0c23da044d8cabe0edd75bc6
  e = 65537

  n = 217037022404812623941372967932339937383440618260689301167305379648414107915459425300548481540493692917022366581760235442444526925713115319185064155523

  _p = 478605014814210760740449684027744153384159299201236688014890786582618769451
  _q = n // _p

  M_ = find_new_M(M, 512)
  print("[+] M' = {}".format(M_))
  p, q = roca_attack(n, M_, 5, 6)
  assert p * q == n
  assert _p in [p, q]
  print("[+] p, q = {}, {}".format(p, q))

if __name__ == "__main__":
  main()

'''
Mon Apr 13 18:43:15 JST 2020 ~/prog/repo/crypto_misc/small_root 96.2697% Full
> time sage roca_attack.sage
[+] M' = 2373273553037774377596381010280540868262890
[+] c' = 1154247
[+] Iteration range: [577123, 1177723]
[+] p, q = 478605014814210760740449684027744153384159299201236688014890786582618769451, 453478370863004904962541767187029241910032539767288755247034216637518065673

real	27m58.600s
user	26m19.136s
sys	0m22.511s
'''
