In the paper [Transcendence of zeros of modular forms](), we investigated when the zeros of certain families of modular forms have transcendental zeros. This repository contains code to compute the non-transcendental zeros of an arbitrary modular form.
We used this code to verify three different conjectures from the paper for small weights. 

- The first conjecture states that other than `i` and `rho`, all the zeros of `Delta_{k,l}` are transcendental. We verified this conjecture for all `k+l <= 1000`.
- The second conjecture states that all of the polynomials `P_{Delta_{k,l}}` are either constant or irreducible. We verified this conjecture for all `k+l <= 1000`.
- The third conjecture states that other than `i` and `rho`, all the zeros of the Miller basis elements are transcendental. We verified this conjecture for all weights `k <= 3000`.


One can run the code as follows:
```
$ sage compute.sage proj_zeros 1000 > out-pz.txt
$ sage compute.sage proj_irred 1000 > out-pi.txt
$ sage compute.sage miller_basis 3000 > out-mb.txt
```
