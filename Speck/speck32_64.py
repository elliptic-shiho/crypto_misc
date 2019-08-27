# An implementation of Speck32/64

# params
k = 16
alpha = 7
beta = 2
MOD = 2**k
MASK = MOD - 1
ROUNDS = 22


def rol(x, y):
  assert 0 < y < k, "Can't shift by negative negative shifts"
  return ((x << y) & MASK) | (x >> (k - y))


def ror(x, y):
  assert 0 < y < k, "Can't shift by negative negative shifts"
  return rol(x, k - y)


def split(x):
  return (x >> k, x & MASK)


def merge(x, y):
  return (x << k) | y


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
  x, y = split(m)
  for i in range(ROUNDS):
    x, y = round_function(x, y, keys[i])
  return merge(x, y)


def decrypt(c, key):
  keys = expand_key(key)
  x, y = split(c)
  for i in range(ROUNDS - 1, -1, -1):
    x, y = round_function_inverse(x, y, keys[i])
  return merge(x, y)


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
  plaintext = 0x6574694c
  ciphertext = 0xa86842f2
  keytext = 0x1918111009080100
  assert encrypt(plaintext, keytext) == ciphertext
  assert decrypt(ciphertext, keytext) == plaintext


if __name__ == "__main__":
  main()
