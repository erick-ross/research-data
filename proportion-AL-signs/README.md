# Proportion of Atkin-Lehner signs

In Theorem 1.1 of [Proportion of Atkin-Lehner sign patterns and Hecke Eigenvalue Equidistribution](), we computed the proportion of modular forms with Atkin-Lehner sign `+1` to be: `1/2 + (1/2) 1_{N = cubefree square} \prod_{p | N} -1/(p^2-p-1) `. In particular, this means that the sign will be biased towards `(-1)^{omega(N)}` whenever `N` is a cubefree square.  In this repository, we provide some data to see this bias numerically.

- `compute-proportions.sage` contains code to compute proportions numerically.
- `proportions.txt` contains the data (obtained by running `sage compute.sage > proportions.txt`).

