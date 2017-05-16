from collections import namedtuple

ECPoint = namedtuple('ECPoint', ['x', 'y', 'z'])

'''
A simple implementation of Elliptic Curve over any ring
'''

class EC(object):
  def __init__(s, F, A, B):
    s.F = F
    s.A = A
    s.B = B

  def is_on_curve(s, P):
    assert isinstance(P, ECPoint)
    x, y, z = s.F(P.x), s.F(P.y), s.F(P.z)
    x3 = s.F.mul(s.F.mul(x, x), x)
    axz2 = s.F.mul(s.F.mul(s.F.mul(s.F(s.A), x), z), z)
    bz3 = s.F.mul(s.F.mul(s.F.mul(s.F(s.B), z), z), z)
    return s.F.add(s.F.add(x3, axz2), bz3) == s.F.mul(s.F.mul(y, y), z)

  def is_infinity(s, P):
    assert isinstance(P, ECPoint)
    return P.x == s.F(0) and P.y == s.F(1) and P.z == s.F(0)

  def add(s, P, Q):
    assert isinstance(P, ECPoint)
    assert isinstance(Q, ECPoint)
    if s.is_infinity(P):
      return Q
    elif s.is_infinity(Q):
      return P
    elif s.F.add(P.y, Q.y) == s.F(0):
      return ECPoint(0, 1, 0)
    Px, Py, Pz = s.F(P.x), s.F(P.y), s.F(P.z)
    Qx, Qy, Qz = s.F(Q.x), s.F(Q.y), s.F(Q.z)
    if s.F.mul(Px, Qy) == s.F.mul(Py, Qx):
      X, Y, Z = Px, Py, Pz
      u = s.F.add(s.F.mul(s.F.mul(s.F(3), X), X), s.F.mul(s.F.mul(s.F(s.A), Z), Z))
      v = s.F.mul(Y, Z)
      a = s.F.mul(Y, v)
      w = s.F.sub(s.F.mul(u, u), s.F.mul(s.F.mul(s.F(8), X), a))
      Rx = s.F.mul(s.F.mul(s.F(2), v), w)
      Ry = s.F.sub(s.F.mul(u, s.F.sub(s.F.mul(s.F.mul(s.F(4), X), a), w)), s.F.mul(s.F.mul(s.F(8), a), a))
      Rz = s.F.mul(s.F.mul(s.F.mul(s.F(8), v), v), v)
    else:
      u = s.F.sub(s.F.mul(Qy, Pz), s.F.mul(Py, Qz))
      v = s.F.sub(s.F.mul(Qx, Pz), s.F.mul(Px, Qz))
      v2 = s.F.mul(v, v)
      v3 = s.F.mul(v2, v)
      w = s.F.sub(s.F.sub(s.F.mul(s.F.mul(s.F.mul(u, u), Pz), Qz), v3), s.F.mul(s.F.mul(s.F.mul(s.F(2), v2), Px), Qz))
      Rx = s.F.mul(v, w)
      Ry = s.F.sub(s.F.mul(u, s.F.sub(s.F.mul(s.F.mul(v2, Px), Qz), w)), s.F.mul(s.F.mul(v3, Py), Qz))
      Rz = s.F.mul(s.F.mul(v3, Pz), Qz)
    return ECPoint(s.F.div(Rx, Rz), s.F.div(Ry, Rz), 1)

  def mul(s, m, P):
    '''
    Scalar Multiplication of P using Binary Method
    '''
    if m < 0:
      b = -1
      m = -m
    else:
      b = 1
    if m == 0:
      return ECPoint(0, 1, 0)
    bits = map(int, bin(m)[2:])[::-1]
    x = P
    if bits[0]:
      res = P
    else:
      res = ECPoint(0, 1, 0)
    for cur in bits[1:]:
      x = s.add(x, x)
      if cur:
        res = s.add(res, x)
    if b == -1:
      res = ECPoint(res.x, s.F.sub(F(0), res.y), res.z)
    return res


