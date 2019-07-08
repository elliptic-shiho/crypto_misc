from vector import Vector
from math import gcd

class IntegerLattice:
  def __init__(s, *args):
    if len(args) == 1 and hasattr(args[0], '__iter__'):
      s.basis = list(args[0])
    else:
      s.basis = args

    if not all(isinstance(v, Vector) for v in s.basis):
      raise ValueError("A lattice basis must be a list of instance of Vector.")
    if not all(len(v) == len(s.basis[0]) for v in s.basis):
      raise ValueError("All lattice basis must have the same size.")
    if not all(all(isinstance(x, int) for x in v) for v in s.basis):
      raise ValueError("This class is only implemented a lattice over the Integer ring.")

    # Initialize "gcd" vector
    v = list(s.basis[0])
    for x in s.basis[1:]:
      x = list(x)
      v = [gcd(s, t) for s, t in zip(x, v)]
    s.gcd_vector = v

  def __repr__(s):
    ret = s.__class__.__name__
    ret += '(' + ', '.join(map(repr, s.basis)) + ')'
    return ret

  def __str__(s):
    return 'Integer Lattice with {} basis [{}]'.format(len(s.basis), ', '.join(map(str, s.basis)))

  def is_point(s, v):
    return all(divmod(x, y)[1] == 0 for x, y in zip(v, s.gcd_vector))

def main():
  bs = [Vector(1, 122, 133, 58, 203)]
  bs += [Vector(0, 259, 0, 0, 0)]
  bs += [Vector(0, 0, 259, 0, 0)]
  bs += [Vector(0, 0, 0, 259, 0)]
  bs += [Vector(0, 0, 0, 0, 259)]

  L = IntegerLattice(bs)

  X = Vector(-4, 30, -14, 27, -35)
  assert L.is_point(X)

  bs = [Vector(1, 0, 0, 0, 12345)]
  bs += [Vector(0, 1, 0, 0, 13333)]
  bs += [Vector(0, 0, 1, 0, 10058)]
  bs += [Vector(0, 0, 0, 1, 1033)]
  bs += [Vector(0, 0, 0, 0, 15432)]

  L = IntegerLattice(bs)

  X = Vector(-2, -3, 5, -1, 0)
  assert L.is_point(X)


main()
