from sage.all import *
import itertools

# display matrix picture with 0 and X
# references: https://github.com/mimoo/RSA-and-LLL-attacks/blob/master/boneh_durfee.sage
def matrix_overview(BB, bound):
  for ii in range(BB.dimensions()[0]):
    a = ('%02d ' % ii)
    for jj in range(BB.dimensions()[1]):
      a += ' ' if BB[ii,jj] == 0 else 'X'
      if BB.dimensions()[0] < 60:
        a += ' '
    if BB[ii, ii] >= bound:
      a += '~'
    print a

def jochemsz_may_trivariate(pol, XX, YY, ZZ, WW, tau, mm):
  '''
  Implementation of Finding roots of trivariate polynomial [1].
  Thanks @Bono_iPad

  References: 
    [1] Ellen Jochemsz and Alexander May. "A Strategy for Finding Roots of Multivariate Polynomials with New Applications in Attacking RSA Variants"
  '''
  tt = floor(mm * tau)
  cond = XX^(7 + 9*tau + 3*tau^2) * (YY*ZZ)^(5+9/2*tau) < WW^(3 + 3*tau)
  print '[+] Bound check: X^{7+9tau+3tau^2} * (YZ)^{5+9/2tau} < W^{3+3tau}:', 
  if cond:
    print 'OK'
  else:
    print 'NG'

  RR = WW * XX^(2*(mm-1)+tt) * (YY*ZZ)^(mm-1)
  # Polynomial constant coefficient (a_0) must be 1
  # XXX: can a_0 be zero?
  f_ = pol
  a0 = f_.constant_coefficient()

  if a0 != 0:
    F = Zmod(RR)
    PK = PolynomialRing(F, 'xs, ys, zs')
    f_ = PR(PK(f_) * F(a0)^-1)

  # Construct set `S` (cf.[1] p.8)
  S = set()
  for i2, i3 in itertools.product(range(0, mm), repeat=2):
    for i1 in range(0, 2*(mm-1) - (i2 + i3) + tt + 1):
      S.add(x^i1 * y^i2 * z^i3)
  S = sorted(S)

  # Construct set `M` (cf.[1] p.8)
  M = set()
  for i2, i3 in itertools.product(range(0, mm + 1), repeat=2):
    for i1 in range(0, 2*mm - (i2 + i3) + tt + 1):
      M.add(x^i1 * y^i2 * z^i3)
  M_S = sorted(M - set(S))
  M = sorted(M)

  # Construct polynomial `g`, `g'` for basis of lattice
  g = []
  g_ = []
  for monomial in S:
    i1 = monomial.degree(x)
    i2 = monomial.degree(y)
    i3 = monomial.degree(z)
    g += [monomial * f_ * XX^(2*(mm-1)+tt-i1) * YY^(mm-1-i2) * ZZ^(mm-1-i3)]

  for monomial in M_S:
    g_ += [monomial * RR]

  # Construct Lattice from `g`, `g'`
  monomials = []
  G = g + g_
  for g_poly in G:
    monomials += g_poly.monomials()
  monomials = sorted(set(monomials))
  assert len(monomials) == len(G)
  dims = len(monomials)
  M = Matrix(IntegerRing(), dims)
  for i in xrange(dims):
    M[i, 0] = G[i](0, 0, 0)
    for j in xrange(dims):
      if monomials[j] in G[i].monomials():
        M[i, j] = G[i].monomial_coefficient(monomials[j]) * monomials[j](XX, YY, ZZ)
  matrix_overview(M, 10)
  print 
  print '=' * 128
  print 

  # LLL

  B = M.LLL()
  matrix_overview(B, 10)

  # Re-construct polynomial `H_i` from Reduced-lattice
  H = [(i, 0) for i in xrange(dims)]
  H = dict(H)
  for j in xrange(dims):
    for i in xrange(dims):
      H[i] += PR((monomials[j] * B[i, j]) / monomials[j](XX, YY, ZZ))

  PX = PolynomialRing(IntegerRing(), 'xn')
  xn = PX.gen()
  PY = PolynomialRing(IntegerRing(), 'yn')
  yn = PX.gen()
  PZ = PolynomialRing(IntegerRing(), 'zn')
  zn = PX.gen()

  # Solve for `x`
  r1 = H[1].resultant(pol, y)
  r2 = H[2].resultant(pol, y)
  r3 = r1.resultant(r2, z)
  x_roots = map(lambda t: t[0], r3.subs(x=xn).roots())
  assert len(x_roots) > 0
  if len(x_roots) == 1 and x_roots[0] == 0:
    print '[-] Can\'t find non-trivial solution for `x`'
    return 0, 0, 0
  x_root = x_roots[0]
  print '[+] Found x0 = %d' % x_root

  # Solve for `z`
  r1_ = r1.subs(x=x_root)
  r2_ = r2.subs(x=x_root)
  z_roots = map(lambda t: t[0], gcd(r1_, r2_).subs(z=zn).roots())
  assert len(z_roots) > 0
  if len(z_roots) == 1 and z_roots[0] == 0:
    print '[-] Can\'t find non-trivial solution for `z`'
    return 0, 0, 0
  z_root = z_roots[0]
  print '[+] Found z0 = %d' % z_root

  # Solve for `y`
  y_roots = map(lambda t: t[0], H[1].subs(x=x_root, z=z_root).subs(y=yn).roots())
  assert len(y_roots) > 0
  if len(y_roots) == 1 and y_roots[0] == 0:
    print '[-] Can\'t find non-trivial solution for `y`'
    return 0, 0, 0
  y_root = y_roots[0]
  print '[+] Found y0 = %d' % y_root
  assert pol(x_root, y_root, z_root) == 0
  return (x_root, y_root, z_root)


if __name__ == '__main__':
  # Sample Implementation: small secret exponents attack for Common Prime RSA proposed at [1]

  # PlaidCTF 2017: Common
  n = 58088432169707511884530887899705328460043341752051203523596713643978064163990486003388043207708974302843002172264417411586749486956628013574926573715095249317804208372132481594468281365459316852927039797580064466731086129431102303726096683810636192790138689214881951836298576801818592777778548978736174404573637727968467102304994883348510912250960513170991711142297420921866513951251779528983034404550484931629702384165892697556071097115581299228167852394258387798734409094231193145682416393190929132395057494573898497863918989076392413537702062176268786330410635086565380004463707094857932051066439903977555824218219175469965671350998705376046031482320586719742092257095632485346554447599297796154806971819544073874759609887168714861438380508381072200369356855424637087904990126639040356052220595703217836390163784968964557550834138598331720685971193948909145336792118311346853933919997080542485529346554810899612986872128053464963893904969598870571593590391792680177336466285483108769678524923709955832516312592396932275858180384073997491565078915852990631405792518132787289319255490058314797418779789436482247156880548274241414788581503832731423362447015575598719157544932430682335826661205441977262084135022109705809490179219483
  e = 1296863142627338816011174237898978985407338763709131570002553514972553044428195580769854555546152719732990573444787592853632636242677007643888487629752955065888264014132275294946697123213980215879592216819850992977582758776020797257958443371755687921634978634111557036569863976062810460472490135115734285081745257140890782551910267112228842110620084653440977034007333330590567591507504424149155838966172077983891168269103390479149309815192817914537419492632845858394205854892694802909512194639800529863998537677351955149803841057222279150345376221210193853988829004456607875070080035867353641697190133291884195010172677474092272163240020759433657106746641878734177642792083528380007586872967990567342937466301144028345682953828338178443155

  gamma = 0.4
  delta = 0.1604

  PR = PolynomialRing(ZZ, 'x, y, z')
  x, y, z = PR.gens()

  # Maximal value of solution `x0`, `y0`, `z0`
  XX = floor(n^delta)
  YY = floor(n^(delta + 0.5 - gamma))
  ZZ = YY

  # Norm of polynomial as vector representation
  WW = floor(n^(2 + 2*delta - 2*gamma))

  # Some Non-negative real (cf. [1] p.13)
  tau = (1/2 + gamma - 4*delta) / (2*delta)

  # Powering degree
  mm = 2

  # Target polynomial
  pol = e^2 * x^2 + e*x*(y+z-2)-(y+z-1)-(n-1)*y*z
  x0, y0, z0 = jochemsz_may_trivariate(pol, XX, YY, ZZ, WW, tau, mm)

  # `x0` is secret exponents. so, `e` * `x0` equivalent to 1 modulo `\phi(n)`.
  assert (Mod(0xdeadbeefcafebabe, n)^e)^x0 == 0xdeadbeefcafebabe

  print '[+] d = %d' % x0

'''
Fri Apr 28 01:55:22 JST 2017 ~/prog/lab/cryptography/crypto-misc/small_root 100%
> time sage jochemsz_may.sage
[+] Bound check: X^{7+9tau+3tau^2} * (YZ)^{5+9/2tau} < W^{3+3tau}: OK
00 X X X X   X X   X X                                                     ~
01   X     X X X       X X   X X                                           ~
02     X     X   X X       X X   X X                                       ~
03       X     X   X X       X X   X X                                     ~
04           X         X   X X         X X   X X                           
05             X         X   X X         X X   X X                         
06                 X         X   X X         X X   X X                     
07                   X         X   X X         X X   X X                   
08                           X           X   X X         X X   X X         
09                             X           X   X X         X X   X X       
10                                 X           X   X X         X X   X X   
11                                   X           X   X X         X X   X X 
12         X                                                               
13               X                                                         
14                     X                                                   
15                       X                                                 
16                         X                                               
17                               X                                         
18                                     X                                   ~
19                                       X                                 ~
20                                         X                               ~
21                                           X                             ~
22                                             X                           ~
23                                               X                         ~
24                                                 X                       ~
25                                                   X                     ~
26                                                     X                   ~
27                                                       X                 ~
28                                                         X               ~
29                                                           X             ~
30                                                             X           ~
31                                                               X         ~
32                                                                 X       ~
33                                                                   X     ~
34                                                                     X   ~
35                                                                       X ~

================================================================================================================================

00 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
01 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
02 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
03 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
04 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
05 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
06 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
07 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
08 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
09 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
10 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
11 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
12 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
13 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
14 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
15 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
16 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
17 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
18 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
19 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
20 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
21 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
22 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
23 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
24 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
25 X                                                                       
26 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
27 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
28 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
29 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
30 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
31 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
32 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
33 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
34 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X ~
35 X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X 
[+] Found x0 = 121409455146529337728771553573777368940787469126222207922323852774150657152520850827361784777080611691463797151113862462191876838854867544454991485396685249135155739665320235632676056272551935226747
[+] Found z0 = 34477692609095619007551451472403824425295545780615560712193766187812541102652373760795720908308312963945151058866515800744120505321466469963393914403824358816801934706230048521867296328521894078409023919964270182962442241107395852721585828400095222273918604839225700476276414948936367757929353286693854439150968251863043
[+] Found y0 = 12378427212203812054196886643302353390095092316143099560226394172278044827874676184190207345725566450126743896447098312828171382944930537650128980339358813086656314388255458990666067762448768584800185978758633519758766969355054481583896558399879737949016525558985089867183042047592892886086277755807216927712208441308304
[+] d = 121409455146529337728771553573777368940787469126222207922323852774150657152520850827361784777080611691463797151113862462191876838854867544454991485396685249135155739665320235632676056272551935226747

real    0m55.516s
user    0m55.280s
sys     0m0.248s
'''
