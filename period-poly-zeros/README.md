# Period Polynomial Zeros
In the paper [Zeros of Even and Odd Period Polynomials](https://arxiv.org/abs/2408.05670), we showed that
(under certain conditions), all the zeros of the even and odd period polynomials
lie on the circle of symmetry `|X|=1/sqrt(N)`.
This repository contains code to support the paper.


The file `period-polys.sage` contains code to compute the zeros of
even/odd period polynomials. One can run the following commands to
verify all of the facts noted in the paper.
```
$ # The file L-VALUES.txt (needed by period-polys.sage) was created
$ # by the command:  "gp --quiet compute-L-values.gp > L-VALUES.txt"
$
$ # [Remark 1.5 (1),(4)]: Verify that kN=12 is optimal
$ sage period-polys.sage L-VALUES.txt optimality_kN_plus
$ sage period-polys.sage L-VALUES.txt optimality_kN_minus
$
$ # [Remark 1.5 (2),(5)]: Verify that betaN is optimal
$ sage period-polys.sage L-VALUES.txt optimality_betaN_plus
$ sage period-polys.sage L-VALUES.txt optimality_betaN_minus
$
$ # [Remark 1.6 (2),(3)]: Verify the claimed locations of exceptional zeros
$ sage period-polys.sage L-VALUES.txt location_values_plus
$ sage period-polys.sage L-VALUES.txt location_values_minus
$
$ # [Conjecture 2.2]: Verify Conjecture 2.2 for k<=100
$ sage period-polys.sage L-VALUES.txt conjecture
$
$ # [Proof of Theorems 1.1,1.2]: Complete the proofs of Theorem 1.1,1.2 (checking small weights)
$ sage period-polys.sage L-VALUES.txt theorem_plus
$ sage period-polys.sage L-VALUES.txt theorem_minus
```


The file `plot.py` contains code to plot several different relevant functions.
This helps to visualize/illustrate several different properties of the functions.
```
$ # [Lemmas 4.1,4.2]: Plot the argument function (illustrating the discontinuity fixes)
$ python3 plot.py arg plus
$ python3 plot.py arg minus
$
$ # [Lemma 4.3]: Plot the radius function (illustrating symmetry and increasing)
$ python3 plot.py rad plus
$ python3 plot.py rad minus
$
$ # [Propositions 5.2,5.3 (1)]: Plot SNMk (illustrating the upper bound)
$ python3 plot.py SNMk plus
$ python3 plot.py SNMk minus
$
$ # [Propositions 5.2,5.3 (2)]: Plot mSNNk (illustrating the upper bound)
$ python3 plot.py mSNNk plus
$ python3 plot.py mSNNk minus
```
