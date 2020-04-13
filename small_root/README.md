Small Root
===========

* [small_root_univariate.sage](small_root_univariate.sage) - Solve univariate modular polynomial with small root (naive method, my practice code, useless).
* [coron.sage](coron.sage) - Solve bivariate integer polynomial with small root (Implementation of [1]).
* [jochemsz_may.sage](jochemsz_may.sage) - Solve trivariate integer polynomial with small root (Implementation of [2]).
* [coppersmith.sage](coppersmith.sage) - Solve univariate modular equation with small root (A *faithfully* Implementation of [3], and I also referenced [4]).
* [howgrave_graham.sage](howgrave_graham.sage) - Solve univariate modular equation with small root (Howgrave-Graham's method. Implementation of [4]).
* [boneh_durfee.sage](boneh_durfee.sage) - Solve bivariate modular equation with small root (Boneh-Durfee's Heuristic Method. Implementation of [5]).
* [focus_group.sage](focus_group.sage) - An implementation of Focus Group Attack[6] against Boneh-Durfee's 0.284 Attack [5].
* [roca_attack.sage](roca_attack.sage) - An implementation of ROCA Attack [7].

# References

* [1] Jean-Sébastien Coron. 2004. [_Finding small roots of bivariate integer polynomial equations revisited._](http://link.springer.com/chapter/10.1007/978-3-540-24676-3_29)
* [2] Ellen Jochemsz and Alexander May. 2006. [_A Strategy for Finding Roots of Multivariate Polynomials with New Applications in Attacking RSA Variants_](http://link.springer.com/chapter/10.1007%2F11935230_18)
* [3] Don Coppersmith. 1996. [_Finding a Small Root of a Univariate Modular Equation_](http://link.springer.com/chapter/10.1007/3-540-68339-9_14)
* [4] Nicholas Howgrave-Graham. 1997. [_Finding Small Roots of Univariate Modular Equations Revisited_](http://link.springer.com/chapter/10.1007/BFb0024458)
* [5] Dan Boneh and Glenn Durfee. 1999. [_Cryptanalysis of RSA with Private Key d Less than N^0.292_](https://link.springer.com/chapter/10.1007/3-540-48910-X_1)
* [6] Stephen D. Miller, Bhargav Narayanan, and Ramarathnam Venkatesan. 2017. [_Coppersmith's lattices and `focus groups': an attack on small-exponent RSA_](https://eprint.iacr.org/2017/835)
* [7] Matus Nemec, Marek Sys, Patr Svenda, Dusan Klinec, and Vashek Matyas. 2017. [_The Return of Coppersmith's Attack: Practical Factorization of Widely Used RSA Moduli_](https://crocs.fi.muni.cz/_media/public/papers/nemec_roca_ccs17_preprint.pdf)

