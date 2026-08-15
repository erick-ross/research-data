# Asymptotics of Hecke polynomial coefficients on the Atkin-Lehner eigenspaces

In the paper [Asymptotics of Hecke polynomial coefficients on the Atkin-Lehner eigenspaces](), we studied the asymptotic behavior of the Hecke polynomial coefficients on the Atkin-Lehner eigenspaces (i.e. the spaces determined by fixing Atkin_lehner sign patterns).
This repository contains code to verify the conjecture stated in the introduction for small parameters.


To run the computations, one can execute the following. 
```
$ # Check k <= 40, N <= 300, m <= 50  over the fullspace and over the newspace 
$ sage compute_coeffs.sage check_conj_c2 40 300 50 FS
$ sage compute_coeffs.sage check_conj_c2 40 300 50 NS
```

