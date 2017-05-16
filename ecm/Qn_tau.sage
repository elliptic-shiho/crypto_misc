'''
A Implementation of Quotient Ring Z/nZ[x] / (x^2 - tau)

 - for Elliptic Curve operation.

References:
* [1] Masaaki Shirase, 2017, "Condition on composite numbers easily factored with elliptic curve method"
'''

class Qn_tau(object):
  """
  Quotient ring of Z/nZ. using irreducible polynomial X^2 - tau
  """
  def __init__(s, n, tau):
    s.n = n
    s.tau = tau

  def add(s, A, B):
    assert isinstance(A, Qn_tau_element)
    assert isinstance(B, Qn_tau_element)
    return Qn_tau_element((A.x0 + B.x0) % s.n, (A.x1 + B.x1) % s.n)

  def sub(s, A, B):
    assert isinstance(A, Qn_tau_element)
    assert isinstance(B, Qn_tau_element)
    return Qn_tau_element((A.x0 - B.x0) % s.n, (A.x1 - B.x1) % s.n)

  def mul(s, A, B):
    assert isinstance(A, Qn_tau_element)
    assert isinstance(B, Qn_tau_element)
    return Qn_tau_element((A.x0 * B.x0 + A.x1 * B.x1 * s.tau) % s.n, (A.x0 * B.x1 + A.x1 * B.x0) % s.n)

  def is_regular(s, A):
    '''
    Is `A` regular?
    cf. [1] Lemma 16
    '''
    assert isinstance(A, Qn_tau_element)
    return gcd(s.n, A.x0^2 - A.x1^2 * s.tau) == 1

  def inv(s, A):
    assert isinstance(A, Qn_tau_element)
    if not s.is_regular(A):
      raise ZeroDivisionError('%s is not regular' % A)
    u = ZZ(inverse_mod(A.x0^2 - A.x1^2 * s.tau, s.n))
    return Qn_tau_element((u * A.x0) % s.n, (u * -A.x1) % s.n)

  def div(s, A, B):
    '''
    Compute A / B
    '''
    assert isinstance(A, Qn_tau_element)
    assert isinstance(B, Qn_tau_element)
    return s.mul(A, s.inv(B))

  def __call__(s, x0, x1=0):
    if isinstance(x0, Qn_tau_element) and x1 == 0:
      return x0
    return Qn_tau_element(x0, x1)

class Qn_tau_element(object):
  def __init__(s, x0, x1):
    s.x0 = x0
    s.x1 = x1

  def __repr__(s):
    return '%s(%r, %r)' % (s.__class__.__name__, s.x0, s.x1)

  def __str__(s):
    return '%d + %dX' % (s.x0, s.x1)

  def __eq__(s, rhs):
    if isinstance(rhs, Qn_tau_element):
      return rhs.x0 == s.x0 and rhs.x1 == s.x1
    else:
      return s.x1 == 0 and rhs == s.x0

