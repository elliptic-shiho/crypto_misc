from __future__ import print_function, division
import unittest

class Vector:
  def __init__(s, *args, **_kwargs):
    kwargs = {
      'base_class': (int, float),
    }
    kwargs.update(_kwargs)
    s.base = kwargs['base_class']

    if len(args) == 1 and hasattr(args[0], '__iter__'):
      s.v = tuple(args[0])
    elif len(args) != 0:
      s.v = tuple(args)
    else:
      s.v = ["dummy"]
    if not all(isinstance(x, s.base) for x in s.v):
      raise ValueError('Invalid Argument Specified: {!r}'.format(args))

  def norm(s):
    import math
    # Calculate euclidean norm
    return math.sqrt(sum(map(lambda x: x**2, s.v)))

  def _prepare_other_as_vector(s, other):
    if not isinstance(other, s.__class__):
      other = s.__class__(other)
    if len(s) != len(other):
      raise ValueError('Vector\'s dimension doesn\'t match: {!r}'.format(other))
    return other

  def add(s, other):
    other = s._prepare_other_as_vector(other)
    return s.__class__(map(lambda x: x[0] + x[1], zip(s.v, other.v)), base_class=s.base)

  def sub(s, other):
    other = s._prepare_other_as_vector(other)
    return s.__class__(map(lambda x: x[0] - x[1], zip(s.v, other.v)), base_class=s.base)

  def scalar_mult(s, other):
    if not isinstance(other, s.base):
      raise ValueError('{!r} must be an instance of {!r}'.format(other, s.base))
    return s.__class__(map(lambda x: other * x, s.v), base_class=s.base)

  def inner_product(s, other):
    try:
      other = s._prepare_other_as_vector(other)
    except ValueError:
      raise ValueError('Inner product can calculate between vector and vector: {!r}'.format(other))
    return sum(map(lambda x: x[0] * x[1], zip(s.v, other.v)))

  def __len__(s):
    return len(s.v)

  def __getitem__(s, n):
    assert n < len(s), 'The specified index out of bounds: {}'.format(n)
    return s.v[n]

  def __iter__(s):
    return iter(s.v)

  def __eq__(s, other):
    if isinstance(other, Vector):
      return s.v == other.v
    elif hasattr(other, '__iter__'):
      return s.v == tuple(other)
    return False

  def __repr__(s):
    return '{}({!r})'.format(s.__class__.__name__, s.v)

  def __str__(s):
    return '({})'.format(', '.join(map(str, s.v)))
