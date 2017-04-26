from sage.all import *
import itertools
import sys
import os
load('AES.sage')

def collect(cipher):
  k = [0] * 15
  res = []
  for i in xrange(256):
    res += [cipher.encrypt([i] + k)]
  return res

def integral_attack(index, cipher, collected_texts):
  cand = []
  for k0 in range(256):
    c = ntopoly(0)
    for i in xrange(256):
      c = c ^^ SBoxInv(collected_texts[i][index] ^^ k0)
    if c == 0:
      cand += [k0]
  return cand

if __name__ == '__main__':
  key = map(ord, os.urandom(16))
  print '[+] Key = %r' % key
  cipher = AES128(4, key, schedule=False)
  print '[+] collecting ciphertexts...',
  sys.stdout.flush()
  ct = collect(cipher)
  print 'done.'
  print '[+] Integral Attack...'
  cand = []
  for i in xrange(16):
    print '[+] i = %d...' % i,
    sys.stdout.flush()
    cand += [integral_attack(i, cipher, ct)]
    print cand[-1]
  print '[+] Key candidate = %r' % cand
  for k in itertools.product(*cand):
    if list(k) == key:
      print '[+] Success'
      break
  else:
    print '[-] Failed'
