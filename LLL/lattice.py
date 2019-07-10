from fractions import Fraction
from vector import Vector
import unittest

def gcd(a, b):
  while b:
    a, b = b, a%b
  return a

class IntegerLattice:
  def __init__(s, *args):
    if len(args) == 1 and hasattr(args[0], '__iter__'):
      s.basis = list(args[0])
    else:
      s.basis = list(args)

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
      v = [gcd(a, b) for a, b in zip(x, v)]
    s.gcd_vector = v

  def __repr__(s):
    ret = s.__class__.__name__
    ret += '(' + ', '.join(map(repr, s.basis)) + ')'
    return ret

  def __str__(s):
    return 'Integer Lattice with {} basis [{}]'.format(len(s.basis), ', '.join(map(str, s.basis)))

  def is_point(s, v):
    return all(divmod(x, y)[1] == 0 for x, y in zip(v, s.gcd_vector))

def gram_schmidt_orthgonalization(L):
  bc = (Fraction, int)
  basis = [Vector(list(x), base_class=bc) for x in L.basis]
  ret = [basis[0]]
  for j in range(1, len(basis)):
    t = Vector([0 for _ in basis], base_class=bc)
    for i in range(j):
      t = t.add(ret[i].scalar_mult(Fraction(basis[j].inner_product(ret[i]), ret[i].inner_product(ret[i]))))
    ret += [basis[j].sub(t)]
  return ret
