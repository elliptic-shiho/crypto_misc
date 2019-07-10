import unittest

from vector import Vector
from lattice import gcd, IntegerLattice, gram_schmidt_orthgonalization

IS_PY2 = not hasattr('dummy', '__iter__')

class TestLattice(unittest.TestCase):
  def test_gcd(s):
    s.assertEqual(gcd(1, 1), 1, 'gcd(1, 1)')
    s.assertEqual(gcd(5, 5), 5, 'gcd(5, 5)')

  def test_basic(s):
    bs1 = [Vector(1, 122, 133, 58, 203)]
    bs1 += [Vector(0, 259, 0, 0, 0)]
    bs1 += [Vector(0, 0, 259, 0, 0)]
    bs1 += [Vector(0, 0, 0, 259, 0)]
    bs1 += [Vector(0, 0, 0, 0, 259)]

    L1 = IntegerLattice(bs1)

    s.assertEqual(L1.basis, bs1)
    s.assertEqual(IntegerLattice(bs1[0], bs1[1], bs1[2], bs1[3], bs1[4]).basis, bs1)

    with s.assertRaises(ValueError) as cm:
      IntegerLattice("hoge")
    s.assertEqual(cm.exception.args[0], "A lattice basis must be a list of instance of Vector.")

    with s.assertRaises(ValueError) as cm:
      IntegerLattice(Vector(0, 1), Vector(1, 2, 3))
    s.assertEqual(cm.exception.args[0], "All lattice basis must have the same size.")

    with s.assertRaises(ValueError) as cm:
      IntegerLattice(Vector(0.5, 1), Vector(0, 1))
    s.assertEqual(cm.exception.args[0], "This class is only implemented a lattice over the Integer ring.")

    s.assertEqual(repr(L1), 'IntegerLattice(Vector((1, 122, 133, 58, 203)), Vector((0, 259, 0, 0, 0)), Vector((0, 0, 259, 0, 0)), Vector((0, 0, 0, 259, 0)), Vector((0, 0, 0, 0, 259)))', '__repr__')

    s.assertEqual(str(L1), 'Integer Lattice with 5 basis [(1, 122, 133, 58, 203), (0, 259, 0, 0, 0), (0, 0, 259, 0, 0), (0, 0, 0, 259, 0), (0, 0, 0, 0, 259)]', '__str__')



  def test_is_point(s):
    bs1 = [Vector(1, 122, 133, 58, 203)]
    bs1 += [Vector(0, 259, 0, 0, 0)]
    bs1 += [Vector(0, 0, 259, 0, 0)]
    bs1 += [Vector(0, 0, 0, 259, 0)]
    bs1 += [Vector(0, 0, 0, 0, 259)]

    L1 = IntegerLattice(bs1)
    X1 = Vector(-4, 30, -14, 27, -35)
    s.assertTrue(L1.is_point(X1), 'X1 is an element of L1')

    bs2 = [Vector(1, 0, 0, 0, 12345)]
    bs2 += [Vector(0, 1, 0, 0, 13333)]
    bs2 += [Vector(0, 0, 1, 0, 10058)]
    bs2 += [Vector(0, 0, 0, 1, 1033)]
    bs2 += [Vector(0, 0, 0, 0, 15432)]

    L2 = IntegerLattice(bs2)

    X2 = Vector(-2, -3, 5, -1, 0)
    s.assertTrue(L2.is_point(X2), 'X2 is an element of L2')

  def test_gram_schmidt(s):
    bs1 = [Vector(1, 122, 133, 58, 203)]
    bs1 += [Vector(0, 259, 0, 0, 0)]
    bs1 += [Vector(0, 0, 259, 0, 0)]
    bs1 += [Vector(0, 0, 0, 259, 0)]
    bs1 += [Vector(0, 0, 0, 0, 259)]

    L1 = IntegerLattice(bs1)
    gs_basis = gram_schmidt_orthgonalization(L1)

    for v1 in gs_basis:
      for v2 in gs_basis:
        if v1 != v2:
          s.assertEqual(v1.inner_product(v2), 0, '')


