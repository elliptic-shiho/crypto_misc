Misc Algorithm Implementations
=====================

* p-1 method [p_1.py](p_1.py)
* Reed-Solomon Code w/ Euclid Decoder Implementation [reed_solomon.sage](reed_solomon.sage)

In `p-1` method, `p-1` must be *`B`-smooth* (i.e. For all integer `x` which satisfy `x | p-1`, It has `x < B).

In this script, we can specify two parameters:

* `B` : maximum bounds of any factor - `B` of `B-smooth`.
* `k` : maximum bounds of power of prime factor of `x`.
<=> If (prime factored `x`, ) `x = p1^e1 p2^e2 ... pt^et` then `max{e1, e2, ..., et} <= k` .
