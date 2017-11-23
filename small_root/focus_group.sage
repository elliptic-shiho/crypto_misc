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

def focus_group_boneh_durfee_284(_pol, modulo, XX, YY, mm, tt, sigma, tau):
  '''
  An implementation of Focus Group Attack [1] against the Boneh-Durfee's "0.284" Attack [2].

  References:
    [1] Stephen D. Miller, Bhargav Narayanan, and Ramarathnam Venkatesan. 2017. "Coppersmith's lattices and `focus groups': an attack on small-exponent RSA"
    [2] Dan Boneh and Glenn Durfee. 1999. "Cryptanalysis of RSA with Private Key d Less than N^0.292"
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
  gil = {}
  # y-shift
  hjl = {}
  # 0 <= \ell <= m
  for l in xrange(0, mm + 1):
    fkmod = f^l * modulo^(mm-l)
    # max(-1, \sigma - \ell) < i <= m - \ell
    for i in xrange(max(0, sigma - l + 1), mm - l + 1):
      gil[i, l] = x^i * fkmod
    # 1 <= j <= min(t, 1 + (\ell - \tau) / 2)
    for j in xrange(1, min(tt + 1, 1 + (l - tau)/2 + 1)):
      hjl[j, l] = y^j * fkmod

  # set of all monomials
  monomials = reduce(lambda u,v: u.union(set(v.monomials())), gil.values(), set())
  monomials = reduce(lambda u,v: u.union(set(v.monomials())), hjl.values(), monomials)

  M = Matrix(ZZ, len(monomials))

  # Convert polynomial to vector
  row = 0
  for l in xrange(0, mm + 1):
    for i in xrange(max(0, sigma - l + 1), mm - l + 1):
      for col, monomial in enumerate(monomials):
        M[row, col] = gil[i, l].monomial_coefficient(monomial) * monomial(XX, YY)
      row += 1
  for l in xrange(0, mm + 1):
    for j in xrange(1, min(tt + 1, 1 + (l - tau)/2 + 1)):
      for col, monomial in enumerate(monomials):
        M[row, col] = hjl[j, l].monomial_coefficient(monomial) * monomial(XX, YY)
      row += 1

  M = M[:row]
  matrix_overview(M)
  B = M.LLL()
  matrix_overview(B)
  Hi = []

  # Generate Polynomials from LLL-reduced Vectors
  for i in xrange(min(5, B.nrows())):
    # print RR(B[i].norm()) < RR(modulo^mm / sqrt(sum([x != 0 for x in B[i]])))
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

def solve_SIP(e, n, delta=0.292, beta=0.5, mm=3, sigma=2, tau=-1):
  F = Zmod(e)
  PR = PolynomialRing(ZZ, 'x, y')
  x, y = PR.gens()
  A = ZZ((n + 1)/2)
  pol = x * (A + y) - 1
  XX = ceil(e ^ delta)
  YY = ceil(e ^ beta)
  tt = floor((1 - 2 * delta) * mm)
  x0, y0 = focus_group_boneh_durfee_284(pol, e, XX, YY, mm, tt, sigma, tau)
  if x0 is None or y0 is None:
    return None
  assert pol(x0, y0) % e == 0
  d0 = (x0 * (A + y0) - 1) / e
  if (Zmod(n)(2)^e)^d0 == 2:
    return d0
  else:
    return -d0

def main():
  # p = 16158503035655503650357438344334975980222051334857742016065172713762327569433945446598600705761456731844358980460949009747059779575245460547544076193224141560315438683650498045875098875194826053398028819192033784138396109321309878080919047169238085235290822926018152521443787945770532904303776199561965192760957166694834171210342487393282284747428088017663161029038902829665513096354230157075129296432088558362971801859230928678799175576150822952201848806616643615613562842355410104862578550863465661734839271290328348967522998634176499319107762583194718667771801067716614802322659239302476074096777926805529798117247
  # q = 16158503035655503650357438344334975980222051334857742016065172713762327569433945446598600705761456731844358980460949009747059779575245460547544076193224141560315438683650498045875098875194826053398028819192033784138396109321309878080919047169238085235290822926018152521443787945770532904303776199561965192760957166694834171210342487393282284747428088017663161029038902829665513096354230157075129296432088558362971801859230928678799175576150822952201848806616643615613562842355410104862578550863465661734839271290328348967522998634176499319107762583194718667771801067716614802322659239302476074096777926805529798117439
  p = 3275342483750763170836416113765340562643588112609761114547434695798746536505772662113665850268902708021591050748320984215116927258714434174724054953133
  q = 3274627040723602338317230751036268460667466921902981431451540870051807157329841903588175940574499055891631204240474172883400239374471379393571624577657
  d = 3001470771525654711865177134747043741463302871182505379927435326735028048350149451
  n = p * q
  e = inverse_mod(d, (p - 1) * (q - 1))
  print d
  d = solve_SIP(e, n, 0.280, mm=7, sigma=3, tau=-1)
  print '[+] d = %d' % d


if __name__ == '__main__':
  main()

'''
Thu Nov 23 17:11:50 JST 2017 ~/prog/lab/cryptography/crypto-misc/small_root 100%
> time sage focus_group.sage
3001470771525654711865177134747043741463302871182505379927435326735028048350149451
00                                                     X                                                                   
01                                                       X                                                                 
02                                                 X                                                                       
03                                                   X                                                                     
04                                                     X     X     X                                                       
05                                                     X X                             X                                   
06                                                 X     X           X                                                     
07                                                 X X           X                                                         
08                                 X                   X   X X     X                                                 X     
09                                                     X X   X     X                 X X                                   
10     X                                           X   X X           X                 X                                   
11                       X                         X X   X       X   X                                                     
12                                 X X                 X   X X X   X                                     X           X   X 
13                                 X                   X X X X     X               X X X                             X     
14     X                               X           X   X X   X     X X               X X                                   
15     X               X X                         X X X X       X   X                 X                                   
16   X                             X X       X         X   X X X   X           X                         X X         X X X 
17                                 X X                 X X X X X   X               X X X       X         X           X   X 
18     X   X                       X   X           X   X X X X     X X             X X X                             X     
19     X               X X   X         X           X X X X   X   X X X               X X                                   
20   X                             X X       X         X X X X X   X           X   X X X     X X         X X         X X X 
21     X   X X                     X X X           X   X X X X X   X X             X X X       X         X           X   X 
22     X   X           X X X X     X   X           X X X X X X   X X X             X X X                             X     
23   X X   X X X                   X X X     X     X   X X X X X   X X         X   X X X     X X         X X         X X X 
24     X   X X         X X X X   X X X X           X X X X X X X X X X             X X X       X         X           X   X 
25   X X   X X X       X X X X X X X X X     X     X X X X X X X X X X         X   X X X     X X         X X         X X X 
26                                                                                               X                         
27                                                                           X X                 X                         
28                                                                         X X                         X                   
29                                                                           X X                 X       X X X             
30                                               X                         X X                         X   X X             
31                                                                           X X                 X       X X X   X   X X X 
32                                               X                         X X                         X   X X X X     X X 
33                                               X                       X X                       X X       X X X X   X   
34                                 X X       X X                   X         X X                 X       X X X   X   X X X 
35                                 X X   X   X X X                         X X                         X   X X X X     X X 
36                                   X   X X X X X                       X X                       X X       X X X X   X   
37                                 X X       X X                   X         X X   X X X   X X X X       X X X   X   X X X 
38                                 X X   X   X X X                         X X     X X   X X X X       X   X X X X     X X 
39                                   X   X X X X X                       X X     X X     X X X X   X X       X X X X   X   
40     X   X X X                   X X X     X X                   X X X     X X   X X X   X X X X       X X X   X   X X X 
41     X   X X X   X               X X X X   X X X                     X   X X     X X   X X X X       X   X X X X     X X 
42         X X X   X X               X X X X X X X                     X X X     X X     X X X X   X X       X X X X   X   
43     X   X X X X     X X X X X X X X X     X X                 X X X X     X X   X X X   X X X X       X X X   X   X X X 
44 X   X   X X X X X   X X X X X X X X X X   X X X                     X   X X     X X   X X X X       X   X X X X     X X 
45 X     X X X X X X X X   X X X X   X X X X X X X                     X X X     X X     X X X X   X X       X X X X   X   

00 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
01 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
02 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
03 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
04 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
05 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
06 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
07 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
08 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
09 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
10 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
11 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
12 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
13 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
14 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
15 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
16 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
17 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
18 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
19 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
20 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
21 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
22 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
23 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
24 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
25 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
26 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
27 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
28 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
29 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
30 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
31 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
32 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
33 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
34 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
35                                                                                               X                         
36                                                                           X X                 X                         
37 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
38 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
39 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
40                                                                           X X                 X       X X X             
41                                                                         X X X                 X     X                   
42 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
43 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
44                                               X                         X X X                 X     X X X X             
45 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 

[+] `x0` = -2022508197860467021722442987001092803918440688862535869295919236688387387115738622
[+] `y0` = -3274984762237182754576823432400804511655527517256371272999487782925276846917807282850920895421700881956611127494397578549258583316592906784147839765395
[+] d = 3001470771525654711865177134747043741463302871182505379927435326735028048350149451

real    0m23.993s
user    0m23.772s
sys     0m0.196s
'''
