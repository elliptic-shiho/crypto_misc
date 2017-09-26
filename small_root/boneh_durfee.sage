from sage.rings.finite_rings.integer_mod_ring import is_IntegerModRing
from sage.all import *
from copy import copy
import itertools
import sys

# display matrix picture with 0 and X
# references: https://github.com/mimoo/RSA-and-LLL-attacks/blob/master/boneh_durfee.sage
def matrix_overview(BB, bound=None):
  for ii in range(BB.dimensions()[0]):
    a = ('%02d ' % ii)
    for jj in range(BB.dimensions()[1]):
      a += ' ' if BB[ii,jj] == 0 else 'X'
      if BB.dimensions()[0] < 60:
        a += ' '
    if bound is not None and BB[ii,ii] > bound:
      a += '~'
    print a
  print 

def boneh_durfee_bivariate(_pol, modulo, XX, YY, mm, tt):
  '''
  Implementaion of Boneh-Durfee's Solving Bivariate Modular Equation with Small Root [1].

  References:
  [1] Dan Boneh and Glenn Durfee. 1999. "Cryptanalysis of RSA with Private Key d Less than N^0.292"
  '''
  '''
  Lifting `_pol` upto Z[x, y]
  '''
  PR = PolynomialRing(ZZ, 'x, y')
  x, y = PR.gens()
  if is_IntegerModRing(_pol.parent().base_ring()):
    '''
    Case `_pol` belongs to (Z/nZ)[x,y]
      Lifting map `h`:
      h: (Z/nZ)[x,y]  -> Z[x',y']
           f(x, y)   |-> f(x', y')
      This map is not homomorphic but we use homomorphism set feature for simplify, in following codes.
    '''
    h = Hom(_pol.parent(), PR)([x, y], check=False) # lift up map
    f = h(_pol)
  else:
    '''
    Case `_pol` belongs to the ring `R[x,y]`
      (includes Z[x,y] and R is an any ring which has natural map to Z).

      In this case, we use the natural map created by SAGE.
    '''
    f = _pol.parent().hom(PR)(_pol)
  # x-shift
  gik = {}
  # y-shift
  hjk = {}
  for k in xrange(mm + 1):
    fkmod = f^k * modulo^(mm-k)
    for i in xrange(mm - k + 1):
      gik[i, k] = x^i * fkmod
    for j in xrange(1, tt + 1):
      hjk[j, k] = y^j * fkmod

  # set of all monomials
  monomials = reduce(lambda u,v: u.union(set(v.monomials())), gik.values(), set())
  monomials = reduce(lambda u,v: u.union(set(v.monomials())), hjk.values(), monomials)
  degx = max([u.degree(x) for u in monomials])
  degy = max([u.degree(y) for u in monomials])

  # sort monomials
  Mx = []
  My = []
  for i in xrange(degx + 1):
    for j in xrange(degy + 1):
      if j > i:
        break
      mono = x^i * y^j
      if mono in monomials:
        Mx += [mono]
  for j in xrange(1, tt):
    for i in xrange(degx + 1):
      jj = i + j
      mono = x^i * y^jj
      if mono in monomials and mono not in Mx:
        My += [mono]
  for j in xrange(degy + 1):
    for i in range(degx + 1):
      mono = x^i * y^j
      if j <= i:
        continue
      if mono in monomials and mono not in Mx and mono not in My:
        My += [mono]
  monomials = Mx + My

  M = Matrix(ZZ, len(monomials))

  # Convert polynomial to vector
  row = 0
  for i in xrange(mm + 1):
    ii = i
    for k in xrange(mm + 1):
      for col, monomial in enumerate(monomials):
        M[row, col] = gik[ii, k].monomial_coefficient(monomial) * monomial(XX, YY)
      row += 1
      if ii == 0:
        break
      ii -= 1
  for j in xrange(1, tt + 1):
    for k in xrange(mm + 1):
      for col, monomial in enumerate(monomials):
        M[row, col] = hjk[j, k].monomial_coefficient(monomial) * monomial(XX, YY)
      row += 1
      if j == 0:
        break

  matrix_overview(M, modulo^mm)
  # remove unhelpful vectors
  for ii in xrange(M.nrows() - 1, len(gik.keys()) - 1, -1):
    if M[ii, ii] > modulo^mm:
      M = M[:ii].stack(M[ii+1:])
  '''
  Gaussian Elimination
  M = M.change_ring(QQ)
  for i in xrange(((mm+1)*(mm+2))/2):
    pivot = M[i, i]
    for j in xrange(i + 1, M.nrows()):
      if M[j, i] != 0:
        M.add_multiple_of_row(j, i, -(M[j, i]/pivot))
  M = M.change_ring(ZZ)
  '''

  matrix_overview(M)
  B = M.LLL()
  matrix_overview(B)
  Hi = []

  # Generate Polynomials from LLL-reduced Vectors
  for i in xrange(min(5, B.nrows())):
    # if RR(B[i].norm()) < RR(modulo^mm / sqrt(len(_pol.monomials()))):
    t_pol = 0
    for j, monomial in enumerate(monomials):
      t_pol += ZZ(B[i, j] / monomial(XX, YY)) * monomial
    Hi += [t_pol]

  # Solve Equation for x
  PK = PolynomialRing(ZZ, 'xk')
  xk = PK.gen()
  root_x = root_y = None
  for h1, h2 in itertools.combinations(Hi, 2):
    h = h1.resultant(h2, y).subs(x=xk)
    if not hasattr(h, 'roots'):
      continue
    roots = h.roots()
    if len(roots) == 0:
      continue
    if len(roots) == 1:
      if roots[0][0] != 0:
        root_x = roots[0][0]
        break
    else:
      for r in roots:
        if r[0] != 0:
          root_x = r[0]
          break
      if root_x is not None:
        break
  if root_x is None:
    print '[-] Can\'t find solution for `x`...'
    return None, None

  print '[+] `x0` = %d' % root_x

  # Solve Equation for y
  for h in Hi:
    h = h.subs(x=root_x, y=xk)
    if not hasattr(h, 'roots'):
      continue
    roots = h.roots()
    if len(roots) == 0:
      continue
    if len(roots) == 1:
      if roots[0][0] != 0:
        root_y = roots[0][0]
        break
    else:
      for r in roots:
        if r[0] != 0:
          root_y = r[0]
          break
      if root_y is not None:
        break
  if root_y is None:
    print '[-] Can\'t find solution for `y`...'
    return None, None

  print '[+] `y0` = %d' % root_y
  return root_x, root_y

def solve_SIP(e, n, delta=0.292, beta=0.5, mm=3):
  F = Zmod(e)
  PR = PolynomialRing(ZZ, 'x, y')
  x, y = PR.gens()
  A = ZZ((n + 1)/2)
  pol = x * (A + y) - 1
  XX = floor(e ^ delta)
  YY = floor(e ^ beta)
  tt = floor((1 - 2 * delta) * mm)
  x0, y0 = boneh_durfee_bivariate(pol, e, XX, YY, mm, tt)
  if x0 is None or y0 is None:
    return None
  assert pol(x0, y0) % e == 0
  d0 = (x0 * (A + y0) - 1) / e
  if (Zmod(n)(2)^e)^d0 == 2:
    return d0
  else:
    return -d0

def main():
  # ASIS CTF Finals 2015: Bodu
  n = 0x3a6160848fb1734cbd0fa22cef582e849223ac04510d51502556b6476d07397f03df155289c20112e87c6f35361d9eb622ca4a0e52d9cd87bf723526c826b88387d06abc4279e353f12ad8ec62ea73c47321a20b89644889a792a73152bc7014b80a693d2e58b123fa925c356b1eba037a4dcac8d8de809167a6fcc30c5c785
  e = 0x365962e8daba7ba92fc08768a5f73b3854f4c79969d5518a078a034437c4669bdb705be4d8b8babf4fda1a6e715269e87b28eecb0d4e02726a27fb8721863740720f583688e5567eb10729bb0d92b322d719949e40c57198d764f1c633e5e277da3d3281ece2ce2eb4df945be5afc3e78498ed0489b2459059664fe15c88a33
  d = solve_SIP(e, n, mm=10, delta=0.28)
  print '[+] d = %d' % d


if __name__ == '__main__':
  main()

'''
Tue Sep 26 01:14:48 JST 2017 ~/prog/lab/crypto/crypto_misc/small_root 100% Full
> time sage boneh_durfee.sage
00 X
01  X                                                                                                            ~
02 XXX
03    X                                                                                                          ~
04  X XX                                                                                                         ~
05 XXXXXX
06       X                                                                                                       ~
07    X  XX                                                                                                      ~
08  X XX XXX
09 XXXXXXXXXX
10           X                                                                                                   ~
11       X   XX                                                                                                  ~
12    X  XX  XXX                                                                                                 ~
13  X XX XXX XXXX
14 XXXXXXXXXXXXXXX
15                X                                                                                              ~
16           X    XX                                                                                             ~
17       X   XX   XXX                                                                                            ~
18    X  XX  XXX  XXXX
19  X XX XXX XXXX XXXXX
20 XXXXXXXXXXXXXXXXXXXXX
21                      X                                                                                        ~
22                X     XX                                                                                       ~
23           X    XX    XXX                                                                                      ~
24       X   XX   XXX   XXXX                                                                                     ~
25    X  XX  XXX  XXXX  XXXXX
26  X XX XXX XXXX XXXXX XXXXXX
27 XXXXXXXXXXXXXXXXXXXXXXXXXXXX
28                             X                                                                                 ~
29                      X      XX                                                                                ~
30                X     XX     XXX                                                                               ~
31           X    XX    XXX    XXXX                                                                              ~
32       X   XX   XXX   XXXX   XXXXX
33    X  XX  XXX  XXXX  XXXXX  XXXXXX
34  X XX XXX XXXX XXXXX XXXXXX XXXXXXX
35 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
36                                     X                                                                         ~
37                             X       XX                                                                        ~
38                      X      XX      XXX                                                                       ~
39                X     XX     XXX     XXXX                                                                      ~
40           X    XX    XXX    XXXX    XXXXX                                                                     ~
41       X   XX   XXX   XXXX   XXXXX   XXXXXX
42    X  XX  XXX  XXXX  XXXXX  XXXXXX  XXXXXXX
43  X XX XXX XXXX XXXXX XXXXXX XXXXXXX XXXXXXXX
44 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
45                                              X                                                                ~
46                                     X        XX                                                               ~
47                             X       XX       XXX                                                              ~
48                      X      XX      XXX      XXXX                                                             ~
49                X     XX     XXX     XXXX     XXXXX                                                            ~
50           X    XX    XXX    XXXX    XXXXX    XXXXXX                                                           ~
51       X   XX   XXX   XXXX   XXXXX   XXXXXX   XXXXXXX
52    X  XX  XXX  XXXX  XXXXX  XXXXXX  XXXXXXX  XXXXXXXX
53  X XX XXX XXXX XXXXX XXXXXX XXXXXXX XXXXXXXX XXXXXXXXX
54 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
55                                                        X                                                      ~
56                                              X         XX                                                     ~
57                                     X        XX        XXX                                                    ~
58                             X       XX       XXX       XXXX                                                   ~
59                      X      XX      XXX      XXXX      XXXXX                                                  ~
60                X     XX     XXX     XXXX     XXXXX     XXXXXX                                                 ~
61           X    XX    XXX    XXXX    XXXXX    XXXXXX    XXXXXXX
62       X   XX   XXX   XXXX   XXXXX   XXXXXX   XXXXXXX   XXXXXXXX
63    X  XX  XXX  XXXX  XXXXX  XXXXXX  XXXXXXX  XXXXXXXX  XXXXXXXXX
64  X XX XXX XXXX XXXXX XXXXXX XXXXXXX XXXXXXXX XXXXXXXXX XXXXXXXXXX
65 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
66                                                                   X                                           ~
67   X                                                               XX                                          ~
68   X XX                                                            XXX                                         ~
69   X XX XXX                                                        XXXX
70   X XX XXX XXXX                                                   XXXXX
71   X XX XXX XXXX XXXXX                                             XXXXXX
72   X XX XXX XXXX XXXXX XXXXXX                                      XXXXXXX
73   X XX XXX XXXX XXXXX XXXXXX XXXXXXX                              XXXXXXXX
74   X XX XXX XXXX XXXXX XXXXXX XXXXXXX XXXXXXXX                     XXXXXXXXX
75   X XX XXX XXXX XXXXX XXXXXX XXXXXXX XXXXXXXX XXXXXXXXX           XXXXXXXXXX
76   X XX XXX XXXX XXXXX XXXXXX XXXXXXX XXXXXXXX XXXXXXXXX XXXXXXXXXXXXXXXXXXXXX
77                                                                              X                                ~
78                                                                    X         XX                               ~
79      X                                                             XX        XXX                              ~
80      X  XX                                                         XXX       XXXX                             ~
81      X  XX  XXX                                                    XXXX      XXXXX                            ~
82      X  XX  XXX  XXXX                                              XXXXX     XXXXXX
83      X  XX  XXX  XXXX  XXXXX                                       XXXXXX    XXXXXXX
84      X  XX  XXX  XXXX  XXXXX  XXXXXX                               XXXXXXX   XXXXXXXX
85      X  XX  XXX  XXXX  XXXXX  XXXXXX  XXXXXXX                      XXXXXXXX  XXXXXXXXX
86      X  XX  XXX  XXXX  XXXXX  XXXXXX  XXXXXXX  XXXXXXXX            XXXXXXXXX XXXXXXXXXX
87      X  XX  XXX  XXXX  XXXXX  XXXXXX  XXXXXXX  XXXXXXXX  XXXXXXXXX XXXXXXXXXXXXXXXXXXXXX
88                                                                                         X                     ~
89                                                                               X         XX                    ~
90                                                                     X         XX        XXX                   ~
91          X                                                          XX        XXX       XXXX                  ~
92          X   XX                                                     XXX       XXXX      XXXXX                 ~
93          X   XX   XXX                                               XXXX      XXXXX     XXXXXX                ~
94          X   XX   XXX   XXXX                                        XXXXX     XXXXXX    XXXXXXX               ~
95          X   XX   XXX   XXXX   XXXXX                                XXXXXX    XXXXXXX   XXXXXXXX
96          X   XX   XXX   XXXX   XXXXX   XXXXXX                       XXXXXXX   XXXXXXXX  XXXXXXXXX
97          X   XX   XXX   XXXX   XXXXX   XXXXXX   XXXXXXX             XXXXXXXX  XXXXXXXXX XXXXXXXXXX
98          X   XX   XXX   XXXX   XXXXX   XXXXXX   XXXXXXX   XXXXXXXX  XXXXXXXXX XXXXXXXXXXXXXXXXXXXXX
99                                                                                                    X          ~
100                                                                                          X         XX         ~
101                                                                                X         XX        XXX        ~
102                                                                      X         XX        XXX       XXXX       ~
103               X                                                      XX        XXX       XXXX      XXXXX      ~
104               X    XX                                                XXX       XXXX      XXXXX     XXXXXX     ~
105               X    XX    XXX                                         XXXX      XXXXX     XXXXXX    XXXXXXX    ~
106               X    XX    XXX    XXXX                                 XXXXX     XXXXXX    XXXXXXX   XXXXXXXX   ~
107               X    XX    XXX    XXXX    XXXXX                        XXXXXX    XXXXXXX   XXXXXXXX  XXXXXXXXX  ~
108               X    XX    XXX    XXXX    XXXXX    XXXXXX              XXXXXXX   XXXXXXXX  XXXXXXXXX XXXXXXXXXX ~
109               X    XX    XXX    XXXX    XXXXX    XXXXXX    XXXXXXX   XXXXXXXX  XXXXXXXXX XXXXXXXXXXXXXXXXXXXXX

00 X
01  X
02 XXX
03    X
04  X XX
05 XXXXXX
06       X
07    X  XX
08  X XX XXX
09 XXXXXXXXXX
10           X
11       X   XX
12    X  XX  XXX
13  X XX XXX XXXX
14 XXXXXXXXXXXXXXX
15                X
16           X    XX
17       X   XX   XXX
18    X  XX  XXX  XXXX
19  X XX XXX XXXX XXXXX
20 XXXXXXXXXXXXXXXXXXXXX
21                      X
22                X     XX
23           X    XX    XXX
24       X   XX   XXX   XXXX
25    X  XX  XXX  XXXX  XXXXX
26  X XX XXX XXXX XXXXX XXXXXX
27 XXXXXXXXXXXXXXXXXXXXXXXXXXXX
28                             X
29                      X      XX
30                X     XX     XXX
31           X    XX    XXX    XXXX
32       X   XX   XXX   XXXX   XXXXX
33    X  XX  XXX  XXXX  XXXXX  XXXXXX
34  X XX XXX XXXX XXXXX XXXXXX XXXXXXX
35 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
36                                     X
37                             X       XX
38                      X      XX      XXX
39                X     XX     XXX     XXXX
40           X    XX    XXX    XXXX    XXXXX
41       X   XX   XXX   XXXX   XXXXX   XXXXXX
42    X  XX  XXX  XXXX  XXXXX  XXXXXX  XXXXXXX
43  X XX XXX XXXX XXXXX XXXXXX XXXXXXX XXXXXXXX
44 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
45                                              X
46                                     X        XX
47                             X       XX       XXX
48                      X      XX      XXX      XXXX
49                X     XX     XXX     XXXX     XXXXX
50           X    XX    XXX    XXXX    XXXXX    XXXXXX
51       X   XX   XXX   XXXX   XXXXX   XXXXXX   XXXXXXX
52    X  XX  XXX  XXXX  XXXXX  XXXXXX  XXXXXXX  XXXXXXXX
53  X XX XXX XXXX XXXXX XXXXXX XXXXXXX XXXXXXXX XXXXXXXXX
54 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
55                                                        X
56                                              X         XX
57                                     X        XX        XXX
58                             X       XX       XXX       XXXX
59                      X      XX      XXX      XXXX      XXXXX
60                X     XX     XXX     XXXX     XXXXX     XXXXXX
61           X    XX    XXX    XXXX    XXXXX    XXXXXX    XXXXXXX
62       X   XX   XXX   XXXX   XXXXX   XXXXXX   XXXXXXX   XXXXXXXX
63    X  XX  XXX  XXXX  XXXXX  XXXXXX  XXXXXXX  XXXXXXXX  XXXXXXXXX
64  X XX XXX XXXX XXXXX XXXXXX XXXXXXX XXXXXXXX XXXXXXXXX XXXXXXXXXX
65 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
66   X XX XXX                                                        XXXX
67   X XX XXX XXXX                                                   XXXXX
68   X XX XXX XXXX XXXXX                                             XXXXXX
69   X XX XXX XXXX XXXXX XXXXXX                                      XXXXXXX
70   X XX XXX XXXX XXXXX XXXXXX XXXXXXX                              XXXXXXXX
71   X XX XXX XXXX XXXXX XXXXXX XXXXXXX XXXXXXXX                     XXXXXXXXX
72   X XX XXX XXXX XXXXX XXXXXX XXXXXXX XXXXXXXX XXXXXXXXX           XXXXXXXXXX
73   X XX XXX XXXX XXXXX XXXXXX XXXXXXX XXXXXXXX XXXXXXXXX XXXXXXXXXXXXXXXXXXXXX
74      X  XX  XXX  XXXX                                              XXXXX     XXXXXX
75      X  XX  XXX  XXXX  XXXXX                                       XXXXXX    XXXXXXX
76      X  XX  XXX  XXXX  XXXXX  XXXXXX                               XXXXXXX   XXXXXXXX
77      X  XX  XXX  XXXX  XXXXX  XXXXXX  XXXXXXX                      XXXXXXXX  XXXXXXXXX
78      X  XX  XXX  XXXX  XXXXX  XXXXXX  XXXXXXX  XXXXXXXX            XXXXXXXXX XXXXXXXXXX
79      X  XX  XXX  XXXX  XXXXX  XXXXXX  XXXXXXX  XXXXXXXX  XXXXXXXXX XXXXXXXXXXXXXXXXXXXXX
80          X   XX   XXX   XXXX   XXXXX                                XXXXXX    XXXXXXX   XXXXXXXX
81          X   XX   XXX   XXXX   XXXXX   XXXXXX                       XXXXXXX   XXXXXXXX  XXXXXXXXX
82          X   XX   XXX   XXXX   XXXXX   XXXXXX   XXXXXXX             XXXXXXXX  XXXXXXXXX XXXXXXXXXX
83          X   XX   XXX   XXXX   XXXXX   XXXXXX   XXXXXXX   XXXXXXXX  XXXXXXXXX XXXXXXXXXXXXXXXXXXXXX
84               X    XX    XXX    XXXX    XXXXX    XXXXXX    XXXXXXX   XXXXXXXX  XXXXXXXXX XXXXXXXXXXXXXXXXXXXXX

00 X
01 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
02 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
04 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
05 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
06 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
07 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
08 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
09 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
10 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
11 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
12 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
13 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
14 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
15 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
16 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
17 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
18 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
19 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
20 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
21 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
22 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
23 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
24 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
25 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
26 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
27 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
28 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
29 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
30 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
31 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
32 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
33 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
34 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
35 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
36 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
37 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
38 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
39 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
40 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
41 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
42 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
43 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
44 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
45 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
46 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
47 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
48 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
49 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
50 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
51 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
52 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
53 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
54 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
55 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
56 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
57 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
58 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
59 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
60 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
61 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
62 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
63 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
64 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
65 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
66 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
67 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
68 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
69 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
70 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
71 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
72 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
73 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
74 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
75 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
76 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
77 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
78 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
79 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
80 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
81 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
82 XXX
83 XXX
84 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

[+] `x0` = -83327572237259130615326494611648145470514336492417406369316732667478074333358086
[+] `y0` = -3241012271559880294449520608078389641936781550113857388795391115360535165237598179564566200793035743956110640188890443265076057478403039028143408511035290
[+] d = 89508186630638564513494386415865407147609702392949250864642625401059935751367507

real  6m33.698s
user  6m27.603s
sys 0m4.012s
'''
