# -*- coding: utf-8 -*-
from sage.all import *

class ReedSolomonCode(object):
  '''
  Reed-Solomon Code Encoder / Decoder Implementation using Euclid Decoder

  References: 
    [1] Shigeichi HIRASAWA and Masao KASAHARA. 2011. "ユークリッド復号法 (Euclidean Decoding Method)" (Japanese paper).
    [2] Yasuo Sugiyama, Masao Kasahara, Shigeichi Hirasawa, and Toshihiko Namekawa. 1975. "A Method for Solving Key Equation for Decoding Goppa Codes"
  '''
  def __init__(s, F, N, K, b=0, t=None, G=None, debug=False):
    '''
    Initialize encoder/decoder.

    Args:
      F     : base finite field (e.g. GF(2^8))
      N     : size of (codeword + message)
      K     : size of message
      b     : (optional) number of burst errors
      t     : (optional) maximum number of terms of error polynomial
      G     : (optional) generator polynomial
      debug : (optional) is print debug message
    '''
    assert N > K > 0
    if t is None:
      t = (N-K)//2
    s.F = F
    s.PR = PolynomialRing(s.F, 'X')
    s.X = s.PR.gen()
    s.zz = s.F.gen()
    s.N = N
    s.K = K
    s.b = b
    s.t = t
    if G is not None:
      s.G = G
    else:
      s.G = reduce(lambda x,y: x*y, [s.X - s.zz^i for i in xrange(b, 2*s.t + b)], 1)
    s.debug = debug

  def to_poly(s, _list):
    '''
    Convert to polynomial from list

    Args:
      _list : An integer-valued list

    Returns:
      A polynomial corresponding to `_list`
    '''
    return sum(map(lambda x: s.X^x[0] * (s.zz^s.F.int_to_log(x[1]) if x[1] != 0 else 0), enumerate(_list)))

  def to_list(s, poly):
    '''
    Convert to list from polynomial

    Args:
      poly : A polynomial over `s.PR`

    Returns:
      An integer-valued list corresponding to `poly`
    '''
    return map(lambda t: s.F.log_to_int(int(t._log_repr())), poly.coefficients(sparse=False))

  def encode(s, message, is_poly=False):
    '''
    Encode message using Reed-Solomon code

    Args:
      message : An integer-valued list
      is_poly : Is return polynomial

    Returns:
      if is_poly == True:
        Reed-Solomon Code polynomial `C(x)`
      else:
        An integer-valued list corresponding to `C(x)`
    '''
    assert len(message) == s.K
    I = s.to_poly(message)
    K = I * s.X^(2*s.t)
    P = K % s.G
    C = K - P # C(x) = k(x) * G(x) = K(x) - (K(x) mod G(x))
    if is_poly:
      return C
    return s.to_list(C)

  def decode(s, Y, is_poly=False):
    '''
    Decode received polynomial using Euclid decoder [1], [2].

    Args:
      Y       : (1) An integer-valued list corresponding to received polynomial `Y(x)`
                (2) A polynomial corresponding to `Y(x)`
      is_poly : Is return polynomial

    Returns:
    if is_poly == True:
      decoded polynomial `C(x)`
    else:
      An integer-valued list corresponsing to decoded message
    '''
    if s.debug:
      print '[+] Y(x):', Y

    if not isinstance(Y, s.PR.Element):
      Y = s.to_poly(Y)
    # Compute Syndrome
    S = []
    for i in xrange(2*s.t):
      S += [Y(s.zz^(i + s.b))]
    s_poly = sum(map(lambda x: x[1] * s.X^x[0], enumerate(S)))
    if s_poly == 0:
      if s.debug:
        print '[+] S(x) = 0'
      return Y
    if s.debug:
      print '[+] S(x):', s_poly

    # Extended Euclidean Algorithm
    # x, y = r[0], r[-1]
    # u, v = A[0], A[-1]
    x, y = s_poly, s.X^(2*s.t)
    u, v = 1, 0
    h = 0
    while x.degree() >= s.t:
      q = y // x
      u, v = v - q * u, u
      x, y = y - q * x, x
      h += 1

    gamma = u(0)^-1
    sigma = gamma * u
    eta = (-1)^h * gamma * x
    if s.debug:
      print '[+] sigma(x):', sigma
      print '[+] eta(x):', eta

    # Key equation check
    assert (sigma * s_poly - eta) % s.X^(2*s.t) == 0, 'key equation check'
    assert (eta * inverse_mod(sigma, s.X^(2*s.t))) % s.X^(2*s.t) == s_poly % (s.X^(2*s.t)), 'key equation check 2'

    # Find root of error-location polynomial
    '''
    Naive Method

    error_pos = []
    for i in xrange(s.F.order()):
      if sigma(s.zz^-i) == 0:
        error_pos += [i]
    '''
    # Use sage's polynoamial root finding algorithm
    error_pos = map(lambda x: s.F.order() - 1 - s.F.int_to_log(x[0].integer_representation()), sigma.roots())
    if s.debug:
      print '[+] Error Position(s):', error_pos

    # Compute derivative of sigma(x)
    sigma_deriv = 0
    for l in error_pos:
      temp = 1
      for i in error_pos:
        if i == l:
          continue
        temp *= (1 - s.zz^i * s.X)
      sigma_deriv += temp * s.zz^l

    assert sigma_deriv == derivative(sigma), 'derivative of sigma'
    if s.debug:
      print '[+] sigma\'(x):', sigma_deriv

    # Find errror-values from error value polynomial by Forney Method
    error_poly = 0
    for i in error_pos:
      error_poly += s.X^i * (-s.zz^(i * (1 - s.b)) * (eta.substitute(X=s.zz^-i) / ((sigma_deriv + s.X).substitute(X=s.zz^-i) - s.zz^-i)))

    if s.debug:
      print '[+] E(x):', error_poly
      print '[+] Corrected Polynomial:', error_poly + Y

    if is_poly:
      return error_poly + Y
    else:
      return s.to_list(error_poly + Y)[s.K:]

  '''
  def decode_lagrange_interpolate(s, Y):
    S = []
    for i in xrange(s.b, 2*s.t + s.b):
      S += [(s.zz^i, Y(s.zz^i))]
    return Y + s.PR.lagrange_polynomial(S)
  '''

def main():
  F = GF(2^8, 'zz', modulo=x^8+x^4+x+3+x^2+1)
  R = ReedSolomonCode(F, 16, 8)
  C = R.encode(map(lambda x: ord(x), 'hogefuga'), True)

  # Sample Error
  E = R.to_poly([0, 0, 20, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 253])

  # Received Polynomial
  Y = C + E
  M2 = R.decode(Y)

  assert ''.join(map(chr, M2)) == 'hogefuga'

if __name__ == '__main__':
  main()
