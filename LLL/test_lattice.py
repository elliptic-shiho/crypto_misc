import unittest

from vector import Vector
from lattice import gcd, IntegerLattice, gram_schmidt_orthgonalization, LLL, is_LLL_basis

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
    L2 = IntegerLattice(bs1[0], bs1[1], bs1[2], bs1[3], bs1[4])

    s.assertEqual(L1.basis, bs1)
    s.assertEqual(L2.basis, bs1)

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

    s.assertEqual(L1, L2)
    s.assertEqual(L1, bs1)
    s.assertNotEqual(L1, 1)

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
          s.assertEqual(v1.inner_product(v2), 0)

  def test_lll(s):
    bs1 = [Vector(4, 1, 2)]
    bs1 += [Vector(4, 7, 2)]
    bs1 += [Vector(3, 1, 7)]

    bs1_expected =  [Vector(4, 1, 2)]
    bs1_expected += [Vector(-1, 0, 5)]
    bs1_expected += [Vector(0, 6, 0)]

    L1 = IntegerLattice(bs1)
    s.assertFalse(is_LLL_basis(L1))

    L2 = LLL(L1)
    s.assertTrue(is_LLL_basis(L2))
    s.assertFalse(is_LLL_basis(L2, delta=1.5)) # delta * ||b[i]|| > ||b[i+1]||

    L2_expected = IntegerLattice(bs1_expected)
    s.assertEqual(L2, L2_expected)

    bs2 = [Vector(1, 122, 133, 58, 203)]
    bs2 += [Vector(0, 259, 0, 0, 0)]
    bs2 += [Vector(0, 0, 259, 0, 0)]
    bs2 += [Vector(0, 0, 0, 259, 0)]
    bs2 += [Vector(0, 0, 0, 0, 259)]

    bs2_expected = [Vector(-4, 30, -14, 27, -35)]
    bs2_expected += [Vector(-23, 43, 49, -39, -7)]
    bs2_expected += [Vector(-13, -32, 84, 23, -49)]
    bs2_expected += [Vector(45, 51, 28, 20, 70)]
    bs2_expected += [Vector(-70, 7, 14, 84, 35)]

    L3 = IntegerLattice(bs2)
    s.assertFalse(is_LLL_basis(L3))

    L4 = LLL(L3)
    s.assertTrue(is_LLL_basis(L4))

    L4_expected = IntegerLattice(bs2_expected)
    s.assertEqual(L4, L4_expected)
