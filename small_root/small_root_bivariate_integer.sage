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

pq, p, q, d = (1759627297431952263480996659875214563L, 741969782135344843L, 2371561942007732041L, 1430131803373864185856432546527493073L)

# Paramaters

P = PolynomialRing(ZZ, 'x, y', order='lex')
x, y = P.gens()

#XX = 2^(16*4)
#YY = 2^(16*4)
XX = 2^(512)
YY = 2^(512)
kk = 2

assert kk >= 0

# Polynomial

f = (33 * x*y + 11 * x + 3 * y + 1) - 16*pq


###############################################################################
##                            algorithm section                              ##
###############################################################################

p00 = f.constant_coefficient()

if gcd(p00, XX * YY) != 1:
  assert False, '[-] Failed: p00 and XY has common divisor'

x, y = f.parent().gens()

print '[+] Bound Check...',
sys.stdout.flush()
WW = ZZ(max(f(x*XX, y*YY).coefficients()))
uu = WW + ((1 - WW) % abs(p00))
delta = max(f.degree(x), f.degree(y))
omega = (delta + kk + 1)^2
nn = uu * (XX * YY)^kk

if RR(XX*YY) < RR(WW^(2/3*f.degree())):
  print 'OK'
else:
  print 'Failed (maybe not found solution...)'

q = (inverse_mod(p00, nn) * f.change_ring(Zmod(nn))).change_ring(ZZ)

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

monomials = map(lambda u: u[1], sorted(map(lambda t: (str(t), t), list(monomials))))
monomials.sort()

M = Matrix(ZZ, omega, omega)
assert len(monomials) == omega

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
print '==='

B = M.LLL()
matrix_overview(B)

PK = PolynomialRing(ZZ, 'xk, yk')
xk, yk = PK.gens()

PX = PolynomialRing(ZZ, 'xs')
xs = PX.gen()

PY = PolynomialRing(ZZ, 'ys')
ys = PY.gen()

monomials = map(lambda t: PK(t), monomials)
pkf = PK(f)

x_root = y_root = None

for i in xrange(omega):
  hi = reduce(lambda u, v: u + v, map(lambda t: ZZ(t[0] / t[1](XX, YY)) * t[1] , zip(B[i], monomials)), 0)
  pol = pkf.resultant(hi, yk).subs(xk=xs)
  if not isinstance(pol, Integer):
    roots = pol.roots()
    if len(roots) > 0:
      if len(roots) == 1 and roots[0][0] == 0:
          continue
      print '[+] Found non-trivial Solution!'
      x_root = roots[0][0]
      print '[+] x0 = %d' % x_root
else:
  if x_root is None:
    print '[-] solution for x was not found...'

for i in xrange(omega):
  hi = reduce(lambda u, v: u + v, map(lambda t: ZZ(t[0] / t[1](XX, YY)) * t[1] , zip(B[i], monomials)), 0)
  pol = pkf.resultant(hi, xk).subs(yk=ys)
  if not isinstance(pol, Integer):
    roots = pol.roots()
    if len(roots) > 0:
      if len(roots) == 1 and roots[0][0] == 0:
          continue
      print '[+] Found non-trivial Solution!'
      y_root = roots[0][0]
      print (3*y_root + 1)/16 == pq
      print '[+] y0 = %d' % y_root
else:
  if y_root is None:
    print '[-] solution for y was not found...'
