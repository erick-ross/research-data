# Code by Dmitriy Shvydkoy and Erick Ross. 

import sys
from sage.schemes.elliptic_curves.cm import is_HCP



ARG = sys.argv[1]
WEIGHT_UB = int(sys.argv[2])
if ARG == 'run_tests': assert WEIGHT_UB >= 132
QEXP_PRECISION = WEIGHT_UB // 12 + 25

# Preliminaries
def E_qexp(k):
    return eisenstein_series_qexp(k, prec=QEXP_PRECISION, normalization='constant')
E4 = E_qexp(4)
E6 = E_qexp(6)
Delta = delta_qexp(prec=QEXP_PRECISION)
j_invar = j_invariant_qexp(prec=QEXP_PRECISION)

x = PolynomialRing(ZZ, 'x').gen()

# Finds j polynomial of given modular form f (given as its q expansion) of weight k
def find_j_poly(f, k):
    if k % 12 == 0:
        a = 0; b = 0; n = k//12
    if k % 12 == 2:
        a = 2; b = 1; n = k//12 - 1
    if k % 12 == 4:
        a = 1; b = 0; n = k//12
    elif k % 12 == 6:
        a = 0; b = 1; n = k//12
    elif k % 12 == 8:
        a = 2; b = 0; n = k//12
    elif k % 12 == 10:
        a = 1; b = 1; n = k//12
    assert k == 12*n + 4*a + 6*b
    
    ret = x-x # zero polynomial
    wt0_form = f / (Delta^n * E4^a * E6^b)
    # Make sure we have enough precision to compute the j poly (plus some extra)
    assert wt0_form.precision_absolute() > 20 

    while wt0_form != 0:
        min_degree = wt0_form.exponents()[0]
        assert min_degree <= 0
        min_coeff = wt0_form[min_degree]
        ret += min_coeff * x^(-min_degree)
        wt0_form -= min_coeff * j_invar^(-min_degree)

    min_deg_f = f.exponents()[0]
    assert ret.degree() == n - min_deg_f
    return ret



# Every modular form of weight k automatically has (a,b) zeros at (rho,i) by 
# the valence formula (the "trivial zeros"). 
# This function returns the discriminants for the nontrivial algebraic zeros of f.
# To find such zeros, one just needs to check which HCP's divide j_poly.
#
# Note that nontrivial zeros at (rho,i) here are counted with multiplicity (1/3,1/2).
# This is because we are using:  f = Delta^l E4^a E6^b P(j). [P = j_poly]
# So because  j(rho) = 0 (of order 3)   and   j(i) = 1728 (of order 2), we have that:
#     ord_0(P) = 1/3 * (ord_rho(f)-a)   and   ord_1728(P) = 1/2 * (ord_i(f)-b).
def find_alg_zeros(f, k):
    j_poly = find_j_poly(f, k)
    # if j_poly == 0: return []
    ret = []
    for irred_factor,mult in j_poly.factor():
        fctr = (irred_factor * irred_factor.denominator()).change_ring(ZZ)
        discr = is_HCP(fctr)
        if discr:
            print(f'FOUND divisor: ({mult}) D={discr} HCP = {fctr} [dividing j_poly = {j_poly}]')
            for i in range(mult): ret.append(discr)
    ret.sort(reverse=True)
    return ret



def cuspidal_proj(k1, k2):
    return E_qexp(k1) * E_qexp(k2) - E_qexp(k1+k2)



#################################################################################
#################################################################################


# Test the function find_alg_zeros() against known results
def run_tests():
    f1 = E_qexp(40) * E_qexp(36) - E_qexp(76)
    assert find_alg_zeros(f1, 76) == []
    # No algebraic zeros.

    f2 = E_qexp(14)^5
    assert find_alg_zeros(f2, 70) == [-3, -3, -3, -4, -4] 
    # At (rho,i); E14 has (2,1) zeros. E14^5 has (10,5) zeros.
    # 70 = 10 mod 12, so E14^5 has (9,4) non-trivial zeros.
    # Counting with multiplicity (1/3,1/2), this means that the
    # j polynomial should detect (3,2) zeros: discrs [-3,-3,-3,-4,-4].

    f3 = Delta^2 * (j_invar^2 + 191025 * j_invar - 121287375)
    assert find_alg_zeros(f3, 24) == [-15]
    # This is carefully constructed so to have zeros with discriminant -15.
    # Here, x^2 + 191025x - 121287375 is the HCP for discriminant -15.

    f4 = f3 + Delta^2 # just change the constant term in the polynomial by 1
    assert find_alg_zeros(f4, 24) == []
    # No algebraic zeros.
    

    MB = victor_miller_basis(126, prec=QEXP_PRECISION)
    f5 = MB[9]
    assert find_alg_zeros(f5, 126) == [-3]
    # g_{126,9} has a non-trivial zero at rho.
    # This is the only non-trivial zero of a Miller basis elt that I could find.
    # Interestingly, g_{126,9} still has all its zeros on the arc (satisfying the conditions of Theorem 1.1).
    

    MB = victor_miller_basis(132, prec=QEXP_PRECISION)
    f6 = MB[9]
    assert find_alg_zeros(f6, 132) == []
    # No algebraic zeros.
    # g_{132,9} is the smallest weight MB example to not have have all its zeros on the lower arc.
    # (So Theorem 1.1 does not apply. But g_{132,9} still happens to have no alg zeros.)
    # See page 2 of https://arxiv.org/pdf/2405.01184 and page 4 of https://www.math.ucla.edu/~wdduke/preprints/serre.pdf
    # Or you can also verify it yourself, by running check_miller_basis_zeros_on_arc().

    for k7 in range(4, WEIGHT_UB+1, 2):
        f7 = E_qexp(k7)
        assert find_alg_zeros(f7, k7) == []
        # No algebraic zeros.
    
    
    print(f'finished run_tests() up to {WEIGHT_UB}')





def check_proj_zeros():
    for k1_k2 in range(8, WEIGHT_UB+1, 2):
        if k1_k2 in [8,10,14]: continue
        for k1 in range(4, (k1_k2-4)+1, 2):
            k2 = k1_k2 - k1
            prj = cuspidal_proj(k1, k2)
            assert find_alg_zeros(prj, k1_k2) == []
    print(f'finished check_proj_zeros() up to {WEIGHT_UB}')


def check_proj_irred():
    for k1_k2 in range(8, WEIGHT_UB+1, 2):
        if k1_k2 in [8,10,14]: continue
        for k1 in range(4, (k1_k2-4)+1, 2):
            k2 = k1_k2 - k1
            prj = cuspidal_proj(k1, k2)
            j_poly = find_j_poly(prj, k1+k2)
            assert j_poly.is_constant() or j_poly.is_irreducible()
    print(f'finished check_proj_irred() up to {WEIGHT_UB}')


def check_miller_basis():
    for k in range(12, WEIGHT_UB+1, 2):
        miller_basis = victor_miller_basis(k, prec=QEXP_PRECISION)
        for m in range(1, len(miller_basis)):
            fm = miller_basis[m]
            alg_zeros = find_alg_zeros(fm, k)
            assert all(abs(D) <= 4 for D in alg_zeros)
            if alg_zeros != []:
                print(f'Weight = {k}, index = {m}: alg_zeros={alg_zeros}' )
    print(f'finished check_miller_basis() up to {WEIGHT_UB}')


    
def check_miller_basis_zeros_on_arc():
    def get_roots(poly):
        real_roots = poly.roots(AA)
        ret = []
        for root, mult in real_roots:
            if 0 <= root <= 1728:
                for i in range(mult):
                    ret.append(root)
        return ret, real_roots

    for k in range(12, WEIGHT_UB+1, 2):
        miller_basis = victor_miller_basis(k, prec=QEXP_PRECISION)
        dim = len(miller_basis) - 1
        for m in range(1, dim+1):
            fm = miller_basis[m]
            j_poly = find_j_poly(fm, k)
            assert j_poly.degree() == dim - m
            good_roots, real_roots = get_roots(j_poly)
            if len(good_roots) != dim - m:
                print(f'FOUND zero not on arc: Weight = {k}, index = {m}, co-index = {dim-m+1}:')
                print(real_roots)
                print(good_roots)



if ARG == 'run_tests':
    run_tests()
elif ARG == 'proj_zeros': 
    check_proj_zeros()
elif ARG == 'proj_irred': 
    check_proj_irred()
elif ARG == 'miller_basis':
    check_miller_basis()
else:
    assert False

