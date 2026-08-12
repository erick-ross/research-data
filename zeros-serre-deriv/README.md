# Interlacing of Zeros of the Serre Derivative of Eisenstein Series

In the paper [Interlacing of Zeros of the Serre Derivative of Eisenstein Series](), we studied interlacing properties of the zeros of Serre derivatives of Eisenstein series.
This repository contains code to check small cases.


To run the computations, one can execute the following. 
```
$ # Verify Stieltjes interlacing for all ell > k
$ sage compute_zeros.sage verify_Stieltjes_interlacing_Fk_Fl
$ # Verify the claimed condition (i.e. Theorem 1.4 and Proposition 5.3) for interlacing.
$ sage compute_zeros.sage verify_interlacing_Fk_Fl
$ # Verify the claimed preliminary estimate (i.e. Lemma 6.1) of zeros of Fsharpk Fharsharpk
$ sage compute_zeros.sage verify_loc_Fsharpk_Fhatsharpk
$ # Verify that the zeros of Fk are interlacing with the zeros of Fsharpk
$ sage compute_zeros.sage check_interlacing_Fk_Fsharpk
$ # Plot Fk - Fhatk and observe that they are within 0.16
$ sage compute_zeros.sage plot_Fk_minus_Fhatk
$ # Plot B0,B1,B2,B3,U,V (and their first,second derivatives)
$  sage compute_zeros.sage plot_B0_B1_B2_B3_U_V
```

