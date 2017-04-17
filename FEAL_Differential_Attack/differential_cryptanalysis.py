from collections import namedtuple, defaultdict
from FEAL import FEAL, XOR, round_function
from multiprocessing import Process, Queue
import itertools
import os

DifferentialPair = namedtuple('DifferentialPair', ['m1', 'm2', 'c1', 'c2'])

def gen_random_bytes(length):
  return map(ord, os.urandom(length))

def make_diff_pair(num, cipher, diff):
  ret = []
  for i in xrange(num):
    m1 = gen_random_bytes(8)
    m2 = XOR(m1, diff)
    c1 = cipher.encrypt(m1)
    c2 = cipher.encrypt(m2)
    ret += [DifferentialPair(m1, m2, c1, c2)]
  return ret

def search_inner(k_cand_base, pairs, delta, queue):
  for k_cand_3 in range(64):
    for k_cand_base2 in itertools.product(range(256), repeat=3):
      score = 0
      k_cand = [16*k_cand_base[0] + k_cand_3, k_cand_base2[0], k_cand_base2[1], k_cand_base2[2]]
      for pair in pairs:
        L1, R1 = pair.c1[:4], pair.c1[4:]
        L2, R2 = pair.c2[:4], pair.c2[4:]
        diff1 = round_function(XOR(R1, k_cand))
        diff2 = round_function(XOR(R2, k_cand))
        outDiff = XOR(XOR(L1, L2), delta)
        diff = XOR(diff1, diff2)
        if diff == outDiff:
          score += 1
        else:
          break
      if score == len(pairs):
        print '[+] Found key: %r' % (k_cand, )
        queue.put(k_cand)
        return

def search_round_delta(pairs, delta):
  print '[+] Searching subkey which has difference %r...' % delta
  q = Queue()
  processes = []
  for k_cand_base1 in range(4):
    p = Process(target=search_inner, args=([k_cand_base1], pairs, delta, q))
    p.start()
  processes += [p]
  k_cand = q.get(True)
  for p in processes:
    p.terminate()
  return k_cand

def search_first(pairs):
  pair = pairs[0]
  R1, L1 = pair.c1[:4], pair.c1[4:]
  R2, L2 = pair.c2[:4], pair.c2[4:]
  mL1, mR1 = pair.m1[:4], pair.m1[4:]
  mL2, mR2 = pair.m2[:4], pair.m2[4:]
  R1 = XOR(L1, R1)
  R2 = XOR(L2, R2)
  k4 = XOR(L1, mL1)
  k5 = XOR(R1, mR1)
  if XOR(L1, mL1) != XOR(L2, mL2) or XOR(R1, mR1) != XOR(R2, mR2):
    return None, None
  return k4, k5


def round_inverse(pairs, k):
  ret = []
  for pair in pairs:
    L1, R1 = pair.c1[:4], pair.c1[4:]
    L2, R2 = pair.c2[:4], pair.c2[4:]
    L1 = XOR(round_function(XOR(R1, k)), L1)
    L2 = XOR(round_function(XOR(R2, k)), L2)
    c1 = R1 + L1
    c2 = R2 + L2
    ret += [DifferentialPair(pair.m1, pair.m2, c1, c2)]
  return ret

def inverse_last(pairs):
  ret = []
  for pair in pairs:
    L1, R1 = pair.c1[:4], pair.c1[4:]
    L2, R2 = pair.c2[:4], pair.c2[4:]
    R1 = XOR(L1, R1)
    R2 = XOR(L2, R2)
    c1 = R1 + L1
    c2 = R2 + L2
    ret += [DifferentialPair(pair.m1, pair.m2, c1, c2)]
  return ret

def main():
  k = [gen_random_bytes(4) for _ in xrange(4 + 2)]
  cipher = FEAL(4, k)
  m = gen_random_bytes(8)
  print '[+] Target Message: %r' % m
  print '[+] Target Subkeys: %r' % (k)
  c = cipher.encrypt(m)
  pair = make_diff_pair(2048, cipher, [0x80, 0x80, 0, 0, 0x80, 0x80, 0, 0])
  pair = inverse_last(pair)
  k3 = search_round_delta(pair, [2, 0, 0, 0])
  pair = round_inverse(pair, k3)
  k2 = search_round_delta(pair, [0x80, 0x80, 0, 0])

  pair = make_diff_pair(2048, cipher, [0x80, 0x80, 0, 0, 0, 0, 0, 0])
  pair = inverse_last(pair)
  pair = round_inverse(pair, k3)
  pair_ = round_inverse(pair, k2)
  k1 = search_round_delta(pair_, [0x80, 0x80, 0, 0])
  if k1 is None:
    pair_ = round_inverse(pair, XOR(k2, [2, 0, 0, 0]))
    k1 = search_round_delta(pair_, [0x80, 0x80, 0, 0])
    k2 = XOR(k2, [2, 0, 0, 0])

  pair = make_diff_pair(2048, cipher, [0x90, 0x90, 0, 0, 0x92, 0x90, 0, 0])
  pair = inverse_last(pair)
  pair = round_inverse(pair, k3)
  pair = round_inverse(pair, k2)
  pair_ = round_inverse(pair, k1)
  k0 = search_round_delta(pair_, [0x90, 0x90, 0, 0])
  if k0 is None:
    pair_ = round_inverse(pair, XOR(k1, [2, 0, 0, 0]))
    k0 = search_round_delta(pair_, [0x90, 0x90, 0, 0])
    k1 = XOR(k1, [2, 0, 0, 0])

  pair = make_diff_pair(2048, cipher, [0, 0, 0, 0, 0, 0, 0, 1])
  pair = inverse_last(pair)
  pair = round_inverse(pair, k3)
  pair = round_inverse(pair, k2)
  pair = round_inverse(pair, k1)
  pair_ = round_inverse(pair, k0)
  k4, k5 = search_first(pair_)
  if k4 is None:
    pair_ = round_inverse(pair, XOR(k0, [2, 0, 0, 0]))
    k4, k5 = search_first(pair_)
    k0 = XOR(k0, [2, 0, 0, 0])
  print '[+] Subkeys: %r' % ([k0, k1, k2, k3, k4, k5])
  cipher_ = FEAL(4, [k0, k1, k2, k3, k4, k5])
  m2 = cipher_.decrypt(c)
  print '[+] Decrypted Ciphertext: %r' % m2
  print '[+] Check: %s' % (m == m2)

if __name__ == '__main__':
  main()
