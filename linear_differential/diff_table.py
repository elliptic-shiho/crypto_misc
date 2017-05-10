from toy_cipher import S
from fractions import gcd

if __name__ == '__main__':
  diffs = [[0 for _ in xrange(16)] for _ in xrange(16)]

  for i in xrange(16):
    for j in xrange(16):
        diffs[i ^ j][S[i] ^ S[j]] += 1

  print '[+] Difference Table: '
  print
  print 'I\O| ' + ' '.join([str(i).ljust(3) for i in xrange(16)])
  print '---+' + '-' * 64
  print '\n'.join(map(lambda x: ('%02x | ' % x[0]) + ' '.join(map(lambda y: str(y).ljust(3), x[1])), enumerate(diffs)))
  print 
  # search maximum difference probability except which has difference 0
  print '[+] Maximum probability:',
  maximum_diff = max([max(t[1:]) for t in diffs[1:]])
  g = gcd(maximum_diff, 16)
  if g == 16:
    print maximum_diff / g 
  else:
    print '%d/%d' % (maximum_diff / g, 16 / g)
