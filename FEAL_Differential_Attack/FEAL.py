def G(x, a, b):
  t = (a + b + x) % 256
  return ((t << 2) & (2**8-1)) | t >> (8-2)

def round_function(x):
  a, b, c, d = x
  b = a ^ b
  c = c ^ d
  b = G(1, b, c)
  c = G(0, c, b)
  a = G(0, a, b)
  d = G(1, d, c)
  return [a, b, c, d]

def XOR(a, b):
  assert len(a) == len(b)
  return map(lambda x: x[0] ^ x[1], zip(a, b))

class FEAL(object):
  def __init__(s, r, keys):
    assert len(keys) == r + 2
    s.r = r
    s.keys = keys

  def encrypt(s, m):
    assert len(m) == 8
    L = XOR(m[:4], s.keys[s.r])
    R = XOR(XOR(m[4:], s.keys[s.r + 1]), L)
    for ROUND in xrange(s.r):
      t = round_function(XOR(R, s.keys[ROUND]))
      L, R = R, XOR(L, t)
    R = XOR(L, R)
    return L + R

  def decrypt(s, m):
    assert len(m) == 8
    L = m[:4]
    R = m[4:]
    R = XOR(L, R)
    L, R = R, L
    for ROUND in xrange(s.r):
      t = round_function(XOR(R, s.keys[s.r - ROUND - 1]))
      L, R = R, XOR(L, t)
    L, R = R, L
    R = XOR(XOR(R, s.keys[s.r + 1]), L)
    L = XOR(L, s.keys[s.r])
    return L + R

if __name__ == '__main__':
  k = [
      [0xde, 0xad, 0xbe, 0xef],
      [0x12, 0x34, 0x56, 0x78],
      [0x87, 0x65, 0x43, 0x21],
      [0xca, 0xfe, 0xba, 0xbe],
      [0x90, 0x90, 0x90, 0x90],
      [0xf0, 0xf0, 0xf0, 0xf0],
  ]
  cipher = FEAL(4, k)
  m = [0xab, 0xcd, 0xef, 0x90, 0x12, 0x34, 0x56, 0x78]
  c = cipher.encrypt(m)
  assert cipher.decrypt(c) == m

