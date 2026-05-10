



def vp(n, p):
    cnt = 0
    while n % p == 0:
        cnt += 1
        n //= p
    return cnt




#######################################3
def psi_local(p,r):
    assert r >= 1
    return (p+1)*p^(r-1)


def psi(n):
    ret = 1
    for p,r in factor(n):
        ret *= psi_local(p,r)
    return ret


################################################




def sigma_local(p, r, alpha):
    assert r >= 1
    assert r >= alpha
    if alpha <= r < 2*alpha:
        return 2 * p^(r-alpha)
    elif r % 2 == 1:
        return 2 * p^((r-1)//2)
    else:
        return (p+1) * p^(r//2-1)


def sigma(n,f):
    ret = 1
    for p,r in factor(n):
        ret *= sigma_local(p,r,vp(f,p))
    return ret


    


#######################################################




def ub_magn_rho_local(p,r,alpha,f_fact):
    assert r >= alpha
    assert r >= 1

    if p == 2:
        return 0
    elif p == 3:
        if r == 1:
            return 1
        else:
            return 0
    elif kronecker(-3,p) == -1:
        return 0
    else:
        return 2
    

# upper bound on the magnitude of rho
def ub_magn_rho(n,f,f_fact):
    ret = 1
    for p,r in factor(n):
        ret *= ub_magn_rho_local(p,r,vp(f,p),f_fact)
    return ret




##################################################3



def ub_magn_rhopm_local(p,r,alpha,f_fact):
    assert r >= alpha
    assert r >= 1

    if p == 2:
        if r == 1:
            return 1
        else:
            return 0
    elif kronecker(-1,p) == -1:
        return 0
    else:
        return 2
        
    

# upper bound on the magnitude of rhopm
def ub_magn_rhopm(n,f,f_fact):
    ret = 1
    for p,r in factor(n):
        ret *= ub_magn_rhopm_local(p,r,vp(f,p),f_fact)
    return ret










def omega(n):
    return len(factor(n))




################################################################
################################################################
#  Here we give code to find all (N,k) that could possibly 
#  have dim S_k(Gamma0(N),chi) <= B  (for some chi)











# For N, f, B given, this function returns kLB where
# dim Sk(Gamma0(N), chi) > B for all k >= kLB  (ind. of chi)
# i.e. we only need to check k < kLB
def find_kLB_f(N, f, B):
    f_fact = factor(f)
    psi_term = psi(N) 
    ub_magn_rho_term = ub_magn_rho(N, f, f_fact)
    ub_magn_rhopm_term = ub_magn_rhopm(N, f, f_fact)
    sigma_term = sigma(N, f)
    magn_const_term = 1
    assert min(psi_term,  ub_magn_rho_term, ub_magn_rhopm_term, sigma_term) >= 0


    for kLB in range(1,100):
        # minimum possible dim over all k >= kLB
        min_dim = (kLB-1)/12 * psi_term 
        min_dim -= (1/3) * ub_magn_rho_term   
        min_dim -= (1/2) * ub_magn_rhopm_term
        min_dim -= (1/2) * sigma_term
        min_dim -= magn_const_term
        if min_dim > B:
            return kLB
    assert False

    

# For N, B given, this function returns kLB where
# dim Sk(Gamma0(N), chi) > B for all k >= kLB  (ind. of chi)
def find_kLB(N, B):
    ret = 1
    for f in divisors(N):
        if not (f % 2 == 0 and N//f % 4 == 2):
            ret = max(ret, find_kLB_f(N, f, B))
    return ret




pairs_N_kLB = []
for N in range(1, 729974):
    kLB = find_kLB(N,2)
    if kLB > 2:
        pairs_N_kLB.append((N,kLB))



print(pairs_N_kLB)


