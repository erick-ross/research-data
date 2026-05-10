# Class Groups of Small Exponent

In the paper, [Distribution of CM points and class groups of small exponent](), we described an algorithm to (conditionally) compute the complete list of (not necessarily fundamental) negative discriminants with class group of small exponent.
This repository contains code to run this algorithm, computing the complete list of discriminants for exponent `1 <= E <= 8`.

The final data is contained in [Table 1](/Table-1.txt), [Table 2](/Table-2.txt), [Table 3](/Table-3.txt), [Table 4](/Table-4.txt), [Table 5](/Table-5.txt), [Table 6](/Table-6.txt), [Table 7](/Table-7.txt), [Table 8](/Table-8.txt).

One can run the code for exponent `E = 8`, for example, by running:

```
$ sage compute-exp.sage 8 > output-8.txt
```

# Computation Notes
See the paper for theoretical justification for details about why the algorithm works.

The code consists of 3 steps:
- `find_fund_discs()`: Compute the complete list of fundamental discriminants with class group of exponent dividing `E`.
- `compute_f_candidates()`: For each fundamental discriminant `D0`, compute the complete list of f candidates (those conductors satisfying the divisibility property from the paper).
- `f_of_expon_E()`: Check all of the f candidates, and return every conductor f with class group of exponent `E`. 

It turns out that only `find_fund_discs()` and `compute_f_candidates() for D0=-3` matter for runtime. These take roughly the following amounts of time to run:

- `E=1,2,3,4` 
  - [< 2 min]    - Everything
- `E=5`  
  - [< 2 min]  - `find_fund_discs()` (0.04M discr bound)
  - [30 min]  - `compute_f_candidates() for D0=-3`   (6M candidates)
- `E=6`
  - [27 min]  - `find_fund_discs()` (6M discr bound)
  - [6 min]    - `compute_f_candidates() for D0=-3`   (1M candidates)
- `E=7`   
  - [< 2 min]  - `find_fund_discs()` (0.1M discr bound)
  - [15 min]  - `compute_f_candidates() for D0=-3`   (3M candidates)
- `E=8`
  - [15 hours] - `find_fund_discs()` (400M discr bound)
  - [3 mins] - `compute_f_candidates() for D0=-3`  (0.8M candidates)







