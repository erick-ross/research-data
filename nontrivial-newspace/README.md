# Nontrivial Newspace

In the paper [Newspaces with Nebentypus: An Explicit Dimension Formula and Classification of Trivial Newspaces](https://arxiv.org/abs/2407.08881), we computed an explicit dimension formula for `S_k^new(Gamma0(N),chi)`. We then used this dimension formula to give lower bounds on `dim S_k^new(Gamma0(N),chi)`  (assuming it is not the case that `2 | f(chi) and 2 || N/f(chi)`; for these `(N,k,chi)`, `dim S_k^new(Gamma0(N),chi) = 0`).
We then used these lower bounds to show that given any bound `B`, `dim S_k^new(Gamma0(N),chi) > B` for `N+k` sufficiently large.
This repository contains code to check the finitely many remaining cases. This allows us to classify all newspaces of small dimension.

The `fullspace` folder contains code to compute the dimension of the full space `S_k(Gamma0(N),chi)` for small `N+k`. The `newspace` folder contains code to compute the dimension of the newspace `S_k^new(Gamma0(N),chi)` for small `N+k`.

In each folder, one can run:
```
$ # verify that our dimension formula matches dimension computed by Sage (for N < 300)
$ sage compute_dim.sage test 300   
$ # compute all of the N, kLB that we need to check
$ sage bound_dim.sage > all_N_kLB.txt   
$ # check all of these N and k < kLB; computing all the triples when dim is small
$ sage compute_dim.sage N_kLB all_N_kLB.txt > ALL_TRIPLES.txt
$ # double check that this gives the same list from just checking all values of N < 200 (and small k).
$ sage compute_dim.sage small_N_k 200 > out_small.txt && diff ALL_TRIPLES.txt out_small.txt
```

- `fullspace/ALL_TRIPLES.txt` gives the complete list of triples `(N,k,chi)` for which `dim S_k(Gamma0(N),chi) <= 2`.
- `newspace/ALL_TRIPLES.txt` gives the complete list of triples `(N,k,chi)` for which `dim S_k^new(Gamma0(N),chi) <= 1` (excluding `(N,k,chi)` where `2 | f(chi) and 2 || N/f(chi)`).
- `disprove_conj.sage` contains code to disprove a conjecture from Greg Martin that `dim S_2^new(Gamma0(N))` takes on all possible non-negative integers.
 
