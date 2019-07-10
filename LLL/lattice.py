from __future__ import division
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

  def __eq__(s, other):
    if isinstance(other, IntegerLattice):
      return s.basis == other.basis
    elif hasattr(other, '__iter__'):
      return s.basis == list(other)
    else:
      return False

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

def is_LLL_basis(L, delta=3/4):
  eps = 1e-10
  n = len(L.basis)
  m = len(L.basis[0])
  gs_basis = gram_schmidt_orthgonalization(L)
  for i in range(1, n):
    bi = L.basis[i]
    for j in range(i):
      bj_star = gs_basis[j]
      if abs(bi.inner_product(bj_star) / (bj_star.norm()**2)) > 1/2:
        return False
  for i in range(n - 1):
    if delta * gs_basis[i].norm()**2 > gs_basis[i+1].norm()**2:
      return False
  return True

def LLL(L, delta=3/4):
  import copy
  import math
  L = copy.deepcopy(L)
  n = len(L.basis)
  while True:
    for i in range(n):
      for _j in range(i):
        bi = L.basis[i]
        j = i - _j - 1
        bj = L.basis[j]
        cij = int(math.floor(bi.inner_product(bj) / bj.inner_product(bj) + 0.5))
        L.basis[i] = bi.sub(bj.scalar_mult(cij))

    gs_basis = gram_schmidt_orthgonalization(L)
    return_flag = True
    for i in range(n - 1):
      if delta * gs_basis[i].norm()**2 > gs_basis[i + 1].norm()**2:
        L.basis[i], L.basis[i + 1] = L.basis[i + 1], L.basis[i]
        return_flag = False
        break
    if return_flag:
      return L

