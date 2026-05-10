# Non-repetition a2 Coefficient

In the paper [Non-repetition of second coefficients of Hecke polynomials](https://arxiv.org/abs/2411.18419), we showed the non-repetition of the second coefficient `a2(Tm(N,2k))` in various aspects.
This repository contains code to compute the a2 coefficient for small `N,k` (among other things). 

- `compute_a2.sage` is code to compute the a2 coefficient.
  -  `sage compute_a2.sage check_decr` verifies that `a2(T2(N,2k))` is decreasing as a function of `k` for `N <= 3392663`. This is used in the proof of Theorem 1.1.
  -  `sage compute_a2.sage check_nonrep_A` and `sage compute_a2.sage check_nonrep_B` verify that `a2(T2(N,2k))` is non-repeating as a function of `k` for `N <= 332427`. This is used in the proof of Theorem 1.2.
  - `sage compute_a2.sage check_T2_T3` verifies that `a2(T3(1,2k)) < a2(T2(1,2k))` for `1 <= k <= 81`. This is used in the proof of Theorem 1.3. 
  - `sage compute_a2.sage check_conj` checks Conjecture 4.1 for `2 <= m <= 10`, `1 <= N <= 100`, `2 <= 2k <= 200`.
- `compute_const.sage` is code to compute the specific constants for the bounds `2^omega(N) <= 4.862  N^(1/4)`  and `sigma0(N) <= 8.447 N^(1/4)`. See [Lemma 2.4 of "Newspaces with Nebentypus: An Explicit Dimension Formula and Classification of Trivial Newspaces"](https://arxiv.org/abs/2407.08881) for a proof of the bounds of this shape. These are used in the proof of Theorem 1.1.
- `compute_one_dim_spaces.sage` is code to compute the eigenvalues of the Hecke operators `T2` and `T4` for one-dimensional spaces in level one (i.e. for `2k = 12,16,18,20,22,26`). 

