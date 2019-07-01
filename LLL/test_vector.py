import unittest

from vector import Vector

IS_PY2 = not hasattr('dummy', '__iter__')

class TestVector(unittest.TestCase):
  def assertAlmostEqual(s, a, b):
    if hasattr(a, '__iter__'):
      a = list(a)
    else:
      a = [a]
    if hasattr(b, '__iter__'):
      b = list(b)
    else:
      b = [b]
    for x in a:
      for y in b:
        super(s.__class__, s).assertAlmostEqual(x, y)

  def test_basic_functions(s):
    # __init__
    v = Vector(1, 2, 3)
    v1 = Vector([1, 2, 3])
    v2 = Vector(0.5, 0.5)

    with s.assertRaises(ValueError) as cm:
      Vector("hoge")
    s.assertEqual(cm.exception.args[0], "Invalid Argument Specified: ('hoge',)")

    with s.assertRaises(ValueError) as cm:
      Vector()
    s.assertEqual(cm.exception.args[0], "Invalid Argument Specified: ()")

    # __str__
    s.assertEqual(str(v), '(1, 2, 3)')
    s.assertEqual(str(v1), '(1, 2, 3)')
    s.assertEqual(str(v2), '(0.5, 0.5)')

    # __repr__
    s.assertEqual(repr(v), 'Vector((1, 2, 3))')

    # __eq__
    s.assertEqual(v, [1, 2, 3])
    s.assertEqual(v, v1)
    s.assertNotEqual(v, 1)
    s.assertEqual([1, 2, 3], v)
    s.assertEqual(v1, v)
    s.assertEqual(tuple(v), (1, 2, 3))

    # __getitem__
    s.assertEqual(v[0], 1)
    s.assertEqual(v[1], 2)
    s.assertEqual(v[2], 3)
    s.assertAlmostEqual(v2, [0.5, 0.5])

    # __iter__
    s.assertEqual(next(iter(v)), 1)

    # __len__
    s.assertEqual(len(v), 3)
    s.assertEqual(len(v2), 2)

  def test_add(s):
    u = Vector(1, 2, 3)
    v = Vector(4, 5, 6)
    w = Vector(0.5, 0.5)

    s.assertEqual(u.add(v), [5, 7, 9])
    s.assertEqual(u.add([4, 5, 6]), [5, 7, 9])

    # consistency check
    s.assertAlmostEqual(w.add(w), [1, 1])
    s.assertAlmostEqual(w, [0.5, 0.5])

    with s.assertRaises(ValueError) as cm:
      u.add([1])
    s.assertEqual(cm.exception.args[0], 'Vector\'s dimension doesn\'t match: Vector((1,))')

  def test_sub(s):
    u = Vector(1, 2, 3)
    v = Vector(4, 5, 6)
    w = Vector(0.5, 0.5)

    s.assertEqual(u.sub(v), [-3, -3, -3])
    s.assertEqual(u.sub([4, 5, 6]), [-3, -3, -3])

    # consistency check
    s.assertAlmostEqual(w.sub(w), [0, 0])
    s.assertAlmostEqual(w, [0.5, 0.5])

  def test_norm(s):
    s.assertAlmostEqual(Vector(1, 2, 3).norm(), 3.7416573867739413)

  def test_scalar_mult(s):
    u = Vector(1, 2, 3)
    v = Vector(4, 5, 6)
    w = Vector(0.5, 0.5)

    s.assertEqual(u.scalar_mult(2), [2, 4, 6])

    # consistency check
    s.assertAlmostEqual(w.scalar_mult(2), [1, 1])

    with s.assertRaises(ValueError) as cm:
      u.scalar_mult("dummy")
    s.assertEqual(cm.exception.args[0], "'dummy' must be an instance of (<{0} 'int'>, <{0} 'float'>)".format(IS_PY2 and "type" or "class"))

  def test_inner_product(s):
    u = Vector(1, 2, 3)
    v = Vector(4, 5, 6)
    w = Vector(0.5, 0.5)

    s.assertEqual(u.inner_product(v), 32)
    s.assertEqual(u.inner_product([4, 5, 6]), 32)

    # consistency check
    s.assertAlmostEqual(w.inner_product(w), [0.5, 0.5])
    s.assertAlmostEqual(w, [0.5, 0.5])

    with s.assertRaises(ValueError) as cm:
      u.inner_product([1])
    s.assertEqual(cm.exception.args[0], "Inner product can calculate between vector and vector: [1]")
