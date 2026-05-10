



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


def beta_psi_f_local(p,r,alpha):
    assert r >= 1
    if alpha == 0:
        if r == 1:
            return p-1
        elif r==2:
            return p^2-p-1
        else:
            return (p^3-p^2-p+1) * p^(r-3)
    else:
        if r==1:
            return p-2
        else:
            return (p^2-2*p+1) * p^(r-2)
    

def beta_psi_f(n, f):
    ret = 1
    for p,r in factor(n):
        ret *= beta_psi_f_local(p,r,vp(f,p))
    return ret


################################################

def magn_beta_sigma_f_local(p,r,alpha):
    assert r>=1
    if alpha == 0:
        if r % 2 == 1:
            return 0
        elif r == 2:
            return p-2
        else:
            return (p^2-2*p+1) * p^(r//2-2)
    elif r == 1:
        if alpha == 1:
            return abs(p-3)/2
        else:
            return p-2
    else:
        if r >= alpha+1 and (r+alpha)%2==1:
            return 0
        elif r >= alpha+2 and (r+alpha)%2==0:
            return (1/2) * (p^2-2*p+1) * p^((r+alpha)/2-2) 
        elif r == alpha:
            return (1/2) * (p^2-3*p+2) * p^(r-2)
        else:
            return (p^2-2*p+1) * p^(r-2) 

# magnitude of beta * sigma_f 
def magn_beta_sigma_f(n,f):
    ret = 1
    for p,r in factor(n):
        ret *= magn_beta_sigma_f_local(p,r,vp(f,p))
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
def ub_magn_rho(n,f_fact):
    vpf = {p:r for (p,r) in f_fact}
    ret = 1
    for p,r in factor(n):
        ret *= ub_magn_rho_local(p,r,vpf[p],f_fact)
    return ret



def magn_beta_rho_f_local(p,r,alpha):
    assert r >= 1
    if r == 1:
        if p == 3 and alpha == 0:
            return 1
        elif p != 3 and alpha >= 1:
            return 1
        elif p != 2 and p != 3 and alpha == 0 and kronecker(-3,p)==1:
            return 0
        else:
            return 2
    elif r == 2:
        if p == 3 and alpha == 0:
            return 1
        elif p != 3 and alpha >= 1:
            return 0
        elif p != 2 and p != 3 and alpha == 0 and kronecker(-3,p)==1:
            return 1
        else:
            return 1
    else:
        if p == 3 and r == 3 and alpha == 0:
            return 1
        else:
            return 0


# magnitude of beta * rho_f
def magn_beta_rho_f(n,f):
    ret = 1
    for p,r in factor(n):
        ret *= magn_beta_rho_f_local(p,r,vp(f,p))
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
        
    

# upper bound on magnitude of rhopm
def ub_magn_rhopm(n,f_fact):
    vpf = {p:r for (p,r) in f_fact}
    ret = 1
    for p,r in factor(n):
        ret *= ub_magn_rhopm_local(p,r,vpf[p],f_fact)
    return ret


def magn_beta_rhopm_f_local(p,r,alpha):
    assert r >= 1
    if r == 1:
        if p == 2 and alpha == 0:
            return 1
        elif p != 2 and alpha >= 1:
            return 1
        elif p != 2 and alpha == 0 and kronecker(-1,p)==1:
            return 0
        else:
            return 2
    elif r == 2:
        if p == 2 and alpha == 0:
            return 1
        elif p != 2 and alpha >= 1:
            return 0
        elif p != 2 and alpha == 0 and kronecker(-1,p)==1:
            return 1
        else:
            return 1
    else:
        if p == 2 and r == 3 and alpha == 0:
            return 1
        else:
            return 0


# magnitude of beta * rhopm_f
def magn_beta_rhopm_f(n,f):
    ret = 1
    for p,r in factor(n):
        ret *= magn_beta_rhopm_f_local(p,r,vp(f,p))
    return ret


def omega(n):
    return len(factor(n))




################################################################
################################################################
#  Here we give code to find all (N,k) that could possibly 
#  have dim S_k^new(Gamma0(N),chi) <= B  (for some chi)











# For N, f, B given, this function returns kLB where
# dim Sk(Gamma0(N), chi) > B for all k >= kLB  (ind. of chi)
def find_kLB_f(N, f, B):
    f_fact = factor(f)
    psi_term = psi(f) * beta_psi_f(N//f, f)
    ub_magn_rho_term = ub_magn_rho(f, f_fact) * magn_beta_rho_f(N//f, f)
    ub_magn_rhopm_term = ub_magn_rhopm(f, f_fact) * magn_beta_rhopm_f(N//f, f)
    magn_sigma_term = 2^omega(f) * magn_beta_sigma_f(N//f, f)
    magn_const_term = abs(moebius(N//f))
    assert min(psi_term,  ub_magn_rho_term, ub_magn_rhopm_term, magn_sigma_term) >= 0


    for kLB in range(1,100):
        # minimum possible dim over all k >= kLB
        min_dim = (kLB-1)/12 * psi_term 
        min_dim -= (1/3) * ub_magn_rho_term   
        min_dim -= (1/2) * ub_magn_rhopm_term
        min_dim -= (1/2) * magn_sigma_term
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
for N in range(1, 424094):
    kLB = find_kLB(N,1)
    if kLB > 2:
        pairs_N_kLB.append((N,kLB))



print(pairs_N_kLB)


