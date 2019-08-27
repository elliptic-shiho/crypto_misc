# An implementation of Speck32/64
from collections import namedtuple, defaultdict
import speck32_64

DiffPair = namedtuple('DifferentialPair', ['p1', 'c1', 'p2', 'c2'])


def main():
  speck32_64.ROUNDS = 3
  plaintext = 0x12345678
  diff = 0x00400000
  keytext = 0x1918111009080100
  epoch = 500
  pairs = []
  round_keys = speck32_64.expand_key(keytext)

  for i in range(epoch):
    plain = (plaintext + i) % (speck32_64.MOD**2)
    c1 = speck32_64.encrypt(plain, keytext)
    c2 = speck32_64.encrypt(plain ^ diff, keytext)
    pairs += [DiffPair(plain, c1, plain ^ diff, c2)]

  cand_dict = defaultdict(int)
  for key_cand in range(256 * 256):
    for pair in pairs:
      c1, c2 = pair.c1, pair.c2
      L1, R1 = speck32_64.split(c1)
      L2, R2 = speck32_64.split(c2)
      L3, R3 = speck32_64.round_function_inverse(L1, R1, key_cand)
      L4, R4 = speck32_64.round_function_inverse(L2, R2, key_cand)
      # if L3 ^ L4 == 0x8000 and R3 ^ R4 == 0x840a:
      if L3 ^ L4 == 0x8100 and R3 ^ R4 == 0x8102:
        cand_dict[key_cand] += 1

  key_3, _ = max(cand_dict.items(), key=lambda x: x[1])
  print(key_3)
  assert key_3 == round_keys[2]


if __name__ == "__main__":
  main()

'''
Tue Aug 27 19:40:25 JST 2019 ~/prog/lab/cryptography/crypto-misc/Speck 100%
> time pypy pure_diff.py
24957

real    0m0.825s
user    0m0.783s
sys     0m0.024s
Tue Aug 27 19:40:28 JST 2019 ~/prog/lab/cryptography/crypto-misc/Speck 100%
> time python pure_diff.py
24957

real    2m6.178s
user    2m6.150s
sys     0m0.056s
'''
