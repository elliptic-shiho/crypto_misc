# An implementation of Speck32/64
from collections import namedtuple, defaultdict

# params
k = 16
alpha = 7
beta = 2
MOD = 2**k
MASK = MOD - 1
ROUNDS = 3

DiffPair = namedtuple('DifferentialPair', ['p1', 'c1', 'p2', 'c2'])


def rol(x, y):
  assert 0 < y < k, "Can't shift by negative negative shifts"
  return ((x << y) & MASK) | (x >> (k - y))


def ror(x, y):
  assert 0 < y < k, "Can't shift by negative negative shifts"
  return rol(x, k - y)


def round_function(x, y, key):
  ret_x = ((ror(x, alpha) + y) % MOD) ^ key
  ret_y = rol(y, beta) ^ ret_x
  return ret_x, ret_y


def round_function_inverse(x, y, key):
  ret_y = ror(x ^ y, beta)
  ret_x = rol(((x ^ key) - ret_y) % MOD, alpha)
  return ret_x, ret_y


def encrypt(m, key):
  keys = expand_key(key)
  # assert len(keys) == ROUNDS, "Invalid keys specified"
  x, y = m >> k, m & MASK
  for i in range(ROUNDS):
    x, y = round_function(x, y, keys[i])
  return (x << k) | y


def decrypt(c, key):
  keys = expand_key(key)
  x, y = c >> k, c & MASK
  for i in range(ROUNDS - 1, -1, -1):
    x, y = round_function_inverse(x, y, keys[i])
  return (x << k) | y


def expand_key(key):
  k_words = []
  while key != 0:
    k_words += [key & MASK]
    key >>= k
  m = len(k_words)
  ret = [k_words[0]]
  ell = k_words[1:]
  for i in range(ROUNDS - 1):
    ell += [((ret[i] + ror(ell[i], alpha)) % MOD) ^ i]
    ret += [rol(ret[i], beta) ^ ell[i + m - 1]]
  return ret


def main():
  plaintext = 0x12345678
  diff = 0x00400000
  keytext = 0x1918111009080100
  epoch = 500
  pairs = []
  round_keys = expand_key(keytext)

  for i in range(epoch):
    plain = (plaintext + i) % (MOD**2)
    c1 = encrypt(plain, keytext)
    c2 = encrypt(plain ^ diff, keytext)
    pairs += [DiffPair(plain, c1, plain ^ diff, c2)]

  cand_dict = defaultdict(int)
  for key_cand in range(256 * 256):
    for pair in pairs:
      c1, c2 = pair.c1, pair.c2
      L1, R1 = c1 >> k, c1 & MASK
      L2, R2 = c2 >> k, c2 & MASK
      L3, R3 = round_function_inverse(L1, R1, key_cand)
      L4, R4 = round_function_inverse(L2, R2, key_cand)
      # if L3 ^ L4 == 0x8000 and R3 ^ R4 == 0x840a:
      if L3 ^ L4 == 0x8100 and R3 ^ R4 == 0x8102:
        cand_dict[key_cand] += 1

  key_3, _ = max(cand_dict.items(), key=lambda x: x[1])
  print(key_3)
  assert key_3 == round_keys[2]


if __name__ == "__main__":
  main()
