from sage.all import *

# Fk := F_2
# F := Fk[x] / (x^8 + x^4 + x^3 + x + 1)
Fk = GF(2)
F.<X> = Fk[x].quotient(x^8+x^4+x^3+x+1)
X = F.gen()


def ntopoly(npoly, X=X):
      return sum(c*X**e for e, c in enumerate(Integer(npoly).bits()))

def polyton(poly, X=X):
  if not hasattr(poly, 'list'):
    if poly in [0, 1]:
      return poly
    poly = poly.polynomial()
  a = poly.list()
  return sum(int(a[i])*(1 << i) for i in xrange(len(a)))


mat_S = Matrix(Fk, 8, 8, [
  [1, 0, 0, 0, 1, 1, 1, 1],
  [1, 1, 0, 0, 0, 1, 1, 1],
  [1, 1, 1, 0, 0, 0, 1, 1],
  [1, 1, 1, 1, 0, 0, 0, 1],
  [1, 1, 1, 1, 1, 0, 0, 0],
  [0, 1, 1, 1, 1, 1, 0, 0],
  [0, 0, 1, 1, 1, 1, 1, 0],
  [0, 0, 0, 1, 1, 1, 1, 1]
])

mat_S_inv = mat_S.inverse()

vec_S = Matrix(Fk, 8, 1, [
  1, 
  1, 
  0, 
  0, 
  0, 
  1, 
  1, 
  0, 
])

vec_S_inv = Matrix(Fk, 8, 1, [
  1, 
  0, 
  1, 
  0, 
  0, 
  0, 
  0, 
  0, 
])

mat_P = Matrix(F, 4, 4, [
  map(ntopoly, [2, 3, 1, 1]),
  map(ntopoly, [1, 2, 3, 1]),
  map(ntopoly, [1, 1, 2, 3]),
  map(ntopoly, [3, 1, 1, 2]),
])

mat_P_inv = mat_P.inverse()

def to_bits(x):
  if hasattr(x, 'list'):
    x = polyton(x)
  return Matrix(Fk, 8, 1, map(ZZ, format(x, '08b'))[::-1])

def to_poly(x):
  return ntopoly(ZZ(''.join(map(lambda t: str(t[0]), x))[::-1], 2))

def SBox(x):
  if not hasattr(x, 'polynomial'):
    x = ntopoly(x)
  if x != 0:
    x = x^-1
  return polyton(to_poly(mat_S * to_bits(x) + vec_S))

def PBox(x):
  x = Matrix(F, 4, 1, map(ntopoly, x))
  return map(lambda t: polyton(t[0]), mat_P * x)

def PBoxInv(x):
  x = Matrix(F, 4, 1, map(ntopoly, x))
  return map(lambda t: polyton(t[0]), mat_P_inv * x)

def SBoxInv(x):
  y = to_poly(mat_S_inv * to_bits(x) + vec_S_inv)
  if y != 0:
    y = y^-1
  return polyton(y)

class AES128(object):
  def __init__(s, rounds, key):
    assert len(key) == 16
    s.key = key
    s.rounds = rounds
    s.keys = []

  def AddRoundKey(s, m, k):
    return m + k

  def MixColumns(s, m):
    return mat_P * m

  def SubBytes(s, m):
    for i in xrange(4):
      for j in xrange(4):
        m[i, j] = ntopoly(SBox(polyton(m[i, j])))
    return m

  def ShiftRows(s, m):
    '''
    [a00, a01, a02, a03]
    [a10, a11, a12, a13]
    [a20, a21, a22, a23]
    [a30, a31, a32, a33]
      =>
    [a00, a01, a02, a03]
    [a11, a12, a13, a10]
    [a22, a23, a20, a21]
    [a33, a30, a31, a32]
    '''
    m = reduce(lambda x, y: x + list(y), Matrix(F, 4, 4, m), [])
    res = m[:4]
    t = m[4:8]
    res += t[1:] + t[:1]
    t = m[8:12]
    res += t[2:] + t[:2]
    t = m[12:16]
    res += t[3:] + t[:3]
    return Matrix(F, 4, 4, res)

  def MixColumnsInv(s, m):
    return mat_P_inv * m

  def SubBytesInv(s, m):
    for i in xrange(4):
      for j in xrange(4):
        m[i, j] = ntopoly(SBoxInv(polyton(m[i, j])))
    return m

  def ShiftRowsInv(s, m):
    '''
    [a00, a01, a02, a03]
    [a10, a11, a12, a13]
    [a20, a21, a22, a23]
    [a30, a31, a32, a33]
      =>
    [a00, a01, a02, a03]
    [a13, a10, a13, a10]
    [a22, a23, a20, a21]
    [a31, a32, a33, a30]
    '''
    m = reduce(lambda x, y: x + list(y), Matrix(F, 4, 4, m), [])
    res = m[:4]
    t = m[4:8]
    res += t[3:] + t[:3]
    t = m[8:12]
    res += t[2:] + t[:2]
    t = m[12:16]
    res += t[1:] + t[:1]
    return Matrix(F, 4, 4, res)

  def KeyExpansions(s):
    keys = []
    # In 128-bit AES, constant value `n` and `b` is 16 and 176.
    # because, 128-bit AES has 10 rounds + 1 AddRoundKey.
    # so, we need 16 * (10 + 1) byte keys to encryption.
    n = 16
    b = (s.rounds + 1) * n
    c = 16

    def rcon(i):
      return polyton(ntopoly(2)^(i-1))

    def rotate(x):
      return x[1:] + x[:1]

    def sch_core(x, i):
      x = rotate(x)
      x = map(polyton, map(SBox, x))
      x[0] = x[0] ^^ rcon(i)
      return x

    key = s.key
    i = 1
    while c < b:
      t = key[-4:]
      if len(key) % 16 == 0:
        t = sch_core(t, i)
        i += 1
      for j in xrange(4):
        key += [key[-16] ^^ t[j]]
        c += 1
      s.keys = key

  def encrypt(s, m):
    assert len(m) == 16

    m = Matrix(F, 4, 4, map(ntopoly, m)).T

    def NextRoundKey():
      t = s.keys[:16]
      s.keys = s.keys[16:]
      return Matrix(F, 4, 4, map(ntopoly, t)).T

    s.KeyExpansions()

    m = s.AddRoundKey(m, NextRoundKey())
    for r in xrange(s.rounds - 1):
      m = s.SubBytes(m)
      m = s.ShiftRows(m)
      m = s.MixColumns(m)
      m = s.AddRoundKey(m, NextRoundKey())
    m = s.SubBytes(m)
    m = s.ShiftRows(m)
    m = s.AddRoundKey(m, NextRoundKey())
    m = reduce(lambda x, y: x + list(y), m.T, [])
    m = map(polyton, m)
    return m

  def decrypt(s, m):
    assert len(m) == 16

    m = Matrix(F, 4, 4, map(ntopoly, m)).T

    def NextRoundKey():
      t = s.keys[-16:]
      s.keys = s.keys[:-16]
      print t, s.keys
      return Matrix(F, 4, 4, map(ntopoly, t)).T

    s.KeyExpansions()

    m = s.AddRoundKey(m, NextRoundKey())
    m = s.ShiftRowsInv(m)
    m = s.SubBytesInv(m)
    for r in xrange(s.rounds - 1):
      m = s.AddRoundKey(m, NextRoundKey())
      m = s.MixColumnsInv(m)
      m = s.ShiftRowsInv(m)
      m = s.SubBytesInv(m)
    m = s.AddRoundKey(m, NextRoundKey())
    m = reduce(lambda x, y: x + list(y), m.T, [])
    m = map(polyton, m)
    return m

if __name__ == '__main__':
  test_vec = [219, 19, 83, 69]
  test_vec_ = [142, 77, 161, 188]
  assert PBox(test_vec) == test_vec_
  assert PBoxInv(test_vec_) == test_vec
  for i in xrange(256):
    assert SBoxInv(SBox(i)) == i

  cipher = AES128(10, [0] * 16)
  cipher.KeyExpansions()
  test_key_expansions = [
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x62, 0x63, 0x63, 0x63, 0x62, 0x63, 0x63, 0x63, 0x62, 0x63, 0x63, 0x63, 0x62, 0x63, 0x63, 0x63,
    0x9b, 0x98, 0x98, 0xc9, 0xf9, 0xfb, 0xfb, 0xaa, 0x9b, 0x98, 0x98, 0xc9, 0xf9, 0xfb, 0xfb, 0xaa,
    0x90, 0x97, 0x34, 0x50, 0x69, 0x6c, 0xcf, 0xfa, 0xf2, 0xf4, 0x57, 0x33, 0x0b, 0x0f, 0xac, 0x99,
    0xee, 0x06, 0xda, 0x7b, 0x87, 0x6a, 0x15, 0x81, 0x75, 0x9e, 0x42, 0xb2, 0x7e, 0x91, 0xee, 0x2b,
    0x7f, 0x2e, 0x2b, 0x88, 0xf8, 0x44, 0x3e, 0x09, 0x8d, 0xda, 0x7c, 0xbb, 0xf3, 0x4b, 0x92, 0x90,
    0xec, 0x61, 0x4b, 0x85, 0x14, 0x25, 0x75, 0x8c, 0x99, 0xff, 0x09, 0x37, 0x6a, 0xb4, 0x9b, 0xa7,
    0x21, 0x75, 0x17, 0x87, 0x35, 0x50, 0x62, 0x0b, 0xac, 0xaf, 0x6b, 0x3c, 0xc6, 0x1b, 0xf0, 0x9b,
    0x0e, 0xf9, 0x03, 0x33, 0x3b, 0xa9, 0x61, 0x38, 0x97, 0x06, 0x0a, 0x04, 0x51, 0x1d, 0xfa, 0x9f,
    0xb1, 0xd4, 0xd8, 0xe2, 0x8a, 0x7d, 0xb9, 0xda, 0x1d, 0x7b, 0xb3, 0xde, 0x4c, 0x66, 0x49, 0x41,
    0xb4, 0xef, 0x5b, 0xcb, 0x3e, 0x92, 0xe2, 0x11, 0x23, 0xe9, 0x51, 0xcf, 0x6f, 0x8f, 0x18, 0x8e]
  assert cipher.keys == test_key_expansions
  cipher = AES128(10, map(ord, '5468617473206D79204B756E67204675'.decode('hex')))
  c = ''.join(map(lambda t: format(t, '02x'), cipher.encrypt(map(ord, '54776F204F6E65204E696E652054776F'.decode('hex')))))
  assert c == '29c3505f571420f6402299b31a02d73a'
  print ''.join(map(lambda t: format(t, '02x'), cipher.decrypt(map(ord, '29c3505f571420f6402299b31a02d73a'.decode('hex')))))
