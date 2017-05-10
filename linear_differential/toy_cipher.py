P = [0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15]
S = [14, 4, 13, 1, 2, 15, 11, 8, 3, 10, 6, 12, 5, 9, 0, 7]
Pinv = [P.index(i) for i in xrange(16)]
Sinv = [S.index(i) for i in xrange(16)]

def split4(x):
  assert 0 <= x < 2**16
  return [x & 0xf, (x >> 4) & 0xf, (x >> 8) & 0xf, (x >> 12) & 0xf]

def merge4(x):
  assert len(x) == 4
  assert all([0 <= t < 16 for t in x])
  return x[0] + x[1] * 2**4 + x[2] * 2**8 + x[3] * 2**12

def split8(x):
  assert 0 <= x < 2**16
  return [x & 0xff, x >> 8]

def merge8(x):
  assert len(x) == 2
  assert all([0 <= t < 256 for t in x])
  return x[0] + x[1] * 2**8

def split16_bits(x):
  assert 0 <= x < 2**16
  return map(int, format(x, '016b')[::-1])

def merge16_bits(x):
  return int(''.join(map(str, x[::-1])), 2)

def Sbox(i):
  return S[i]

def Pbox(i):
  i = split16_bits(i)
  o = split16_bits(0)
  for j in xrange(16):
    o[P[j]] = i[j]
  return merge16_bits(o)

def SboxInv(i):
  return Sinv[i]

def PboxInv(i):
  i = split16_bits(i)
  o = split16_bits(0)
  for j in xrange(16):
    o[Pinv[j]] = i[j]
  return merge16_bits(o)

class ToyCipher(object):
  def __init__(s, key, rounds=4):
    assert 0 <= key < 2**16
    s.k = key
    s.r = rounds

  def encrypt(s, m):
    assert 0 <= m < 2**16
    for r in xrange(s.r):
      m = Pbox(merge4(map(Sbox, split4(m)))) ^ s.k
    return m

  def decrypt(s, c):
    assert 0 <= c < 2**16
    for r in xrange(s.r):
      c = merge4(map(SboxInv, split4(PboxInv(c ^ s.k))))
    return c

def main():
  cipher = ToyCipher(0xdead)
  m = 0x1234
  c = cipher.encrypt(m)
  m2 = cipher.decrypt(c)
  print hex(c)
  print m2 == m

if __name__ == '__main__':
  main()
