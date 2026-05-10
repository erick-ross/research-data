import math
import sys


def get_classnum(D):
    return BQFClassGroup(D).order()

def get_expon(D):
    return BQFClassGroup(D).abelian_group().exponent()

def get_group_str(D):
    invars = BQFClassGroup(D).abelian_group().invariants()
    return ' x '.join(f'Z{d}' for d in invars)



# Find all fundamental discriminants of exponent dividing E (uses known bounds)
def find_fund_discs(E):
    # bounds from https://arxiv.org/pdf/1803.02056
    if   E==1: bound = 163
    elif E==2: bound = 5460
    elif E==3: bound = 4027
    elif E==4: bound = 435435
    elif E==5: bound = 37363
    elif E==6: bound = 5761140
    elif E==7: bound = 118843
    elif E==8: bound = 430950520


    ret = []
    for D0_ in range(-3, -bound-1, -1):
        D0 = ZZ(D0_)
        if not D0.is_fundamental_discriminant(): continue

        # Since sqfr(|G|) | exp(G), we can continue early if sqfr(cn) does not divide E
        cn = get_classnum(D0)
        if E % squarefree_part(cn) != 0: continue 
        
        # Now, actually check if D0 exponent divides E 
        if E % get_expon(D0) == 0:
            ret.append(D0)
    return ret

# evalutate L(p^r) for fund discr D0
def get_L_pr(p, r, D0):
    symb = kronecker(D0, p)
    return p^(r-1) * (p - symb)

# evaluate L(f) for fund discr D0
def get_L(f, D0):
    ret = 1
    for p,r in list(factor(f)):
        ret = lcm(ret, get_L_pr(p, r, D0))
    return ret

# Compute all prime powers p^r <= X
def get_prime_powers(X):
    ret = []
    for p in prime_range(X+1):
        r = 1
        while p^r <= X:
            ret.append((p,r))
            r += 1
    return ret

def theta(m, D0):   
    ret = 1
    for p,r in get_prime_powers(2*m):
        if m % get_L_pr(p, r, D0) == 0:
            ret *= p^r 
    return ret

# Now, we extend to discriminants D = D0*f^2 with conductor f >= 1.
# Let u_{-3} = 6, u_{-4} = 4, and u_D0 = 2 for D0 <= -7
# We know that if O_D has exponent E, then L(f) | 12 E u_D0. (Lemma)
# So we compute the complete list of such "f candidates".

# compute the complete list of f candidates: the f such that L(f) | 12 E u_D0
def compute_f_candidates(D0, E):
    if   D0 == -3: u_D0 = 6
    elif D0 == -4: u_D0 = 4
    else:          u_D0 = 2
    _12_E_uD0 = 12 * E * u_D0
    # We know that  f | theta(L(f)). (Lemma)
    # So the only possible candidates f are divisors of theta(m) for m | 12 E u_D0.
    f_precandidates = set()
    for L in divisors(_12_E_uD0):
        for f in divisors(theta(L, D0)):
            f_precandidates.add(f)

    cands = []
    for f in f_precandidates:
        if _12_E_uD0  % get_L(f, D0) == 0:
            cands.append(f)
    cands.sort()
    return cands


# Return the f in f_vals with exponent E.
# This is highly optimized by caching and early return.
def f_of_expon_E(D0, f_vals, E):
    EXPON_COND = {}
    def expon_cond(f):
        if f not in EXPON_COND:
            EXPON_COND[f] = get_expon(D0*f^2)
        return EXPON_COND[f]
    
    def has_expon_E(f):
        # Check small divisors f' of f first.
        # This will exclude essentially all the f with exponent != E.
        for f_pm in range(2, 1000):
            if f % f_pm != 0: continue
            if E % expon_cond(f_pm) != 0:
                return False
        # Finally, return the correct answer
        return expon_cond(f) == E
    
    return [f for f in f_vals if has_expon_E(f)]
            



# compute the complete list of discriminants with exponent E
def compute_discrs_with_exp(E):
    print(f'Exponent E = {E} ###########')
    fund_discs = find_fund_discs(E)
    print(f'fundamental discriminants with exponent dividing E: {len(fund_discs)}')
    print(fund_discs)

    sol = []
    fund_discs.reverse()
    for D0 in fund_discs:
        f_cands = compute_f_candidates(D0, E)
        print(f'fund_discr = {D0}: {len(f_cands)} f_cands, all <= {max(f_cands)}')
        f_expon_E = f_of_expon_E(D0, f_cands, E)
        for f in f_expon_E:
            sol.append( (D0, f, D0*f^2, get_group_str(D0*f^2)) )
    sol.sort(key=lambda x: (abs(x[0]), x[1]))

    print()
    print(f'All discriminants with exponent E={E}:')
    print(len(sol))
    print(sol)
    lst = []
    dic = {}
    for tup in sol:
        D = tup[2]
        lst.append(D)
        gp_str = tup[3]
        if gp_str not in dic:
            dic[gp_str] = []
        dic[gp_str].append(D)

    print()
    print('Discriminant list:')
    print(f'Length = {len(lst)},   Largest value = {min(lst)}')
    print(sorted(lst, reverse=True))
    print()
    print('Discriminants by Group:')
    for ky in dic:
        print(f'Group={repr(ky)}: {sorted(dic[ky], reverse=True)}')
    print()
    
        





# Testing the possible conductors:
ARG = ZZ(int(sys.argv[1]))
assert 1 <= ARG <= 8
compute_discrs_with_exp(ARG)

    

    