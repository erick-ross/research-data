import sys

ARG = sys.argv[1]
ARG_PARAM = sys.argv[2]




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

def beta_sigma_f_local(p,r,alpha):
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
            return (p-3)/2
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


def beta_sigma_f(n,f):
    ret = 1
    for p,r in factor(n):
        ret *= beta_sigma_f_local(p,r,vp(f,p))
    return ret

#######################################################

def get_chi_p_alpha(f_fact, chi, p, x):
    rems = [(int(x) if q == p else 1) for (q,alpha) in f_fact]
    mods = [p^alpha for (p,alpha) in f_fact]
    assert rems.count(1) == len(rems) - 1
    x_hat = CRT(rems, mods)
    chi_prim = chi.primitive_character()
    return chi_prim(x_hat)




def rho_local(p,r,alpha,chi,f_fact):
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
        u = mod(-3,p^r).sqrt()
        chi_x = get_chi_p_alpha(f_fact, chi, p, (-1+u)/2)
        assert chi_x^3 == 1
        if chi_x == 1:
            return 2
        else:
            return -1
    


def rho(n,chi,f_fact):
    vpf = {p:r for (p,r) in f_fact}
    ret = 1
    for p,r in factor(n):
        ret *= rho_local(p,r,vpf[p],chi,f_fact)
    return ret


def beta_rho_f_local(p,r,alpha):
    assert r >= 1
    if r == 1:
        if p == 3 and alpha == 0:
            return -1
        elif p != 3 and alpha >= 1:
            return -1
        elif p != 2 and p != 3 and alpha == 0 and kronecker(-3,p)==1:
            return 0
        else:
            return -2
    elif r == 2:
        if p == 3 and alpha == 0:
            return -1
        elif p != 3 and alpha >= 1:
            return 0
        elif p != 2 and p != 3 and alpha == 0 and kronecker(-3,p)==1:
            return -1
        else:
            return 1
    else:
        if p == 3 and r == 3 and alpha == 0:
            return 1
        else:
            return 0



def beta_rho_f(n,f):
    ret = 1
    for p,r in factor(n):
        ret *= beta_rho_f_local(p,r,vp(f,p))
    return ret



##################################################3



def rhopm_local(p,r,alpha,chi,f_fact):
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
        upm = mod(-1,p^r).sqrt()
        chi_x = get_chi_p_alpha(f_fact, chi, p, upm)
        assert chi_x^4 == 1
        if chi_x == 1:
            return 2
        elif chi_x == -1:
            return -2
        else:
            return 0
        
    


def rhopm(n,chi,f_fact):
    vpf = {p:r for (p,r) in f_fact}
    ret = 1
    for p,r in factor(n):
        ret *= rhopm_local(p,r,vpf[p],chi,f_fact)
    return ret


def beta_rhopm_f_local(p,r,alpha):
    assert r >= 1
    if r == 1:
        if p == 2 and alpha == 0:
            return -1
        elif p != 2 and alpha >= 1:
            return -1
        elif p != 2 and alpha == 0 and kronecker(-1,p)==1:
            return 0
        else:
            return -2
    elif r == 2:
        if p == 2 and alpha == 0:
            return -1
        elif p != 2 and alpha >= 1:
            return 0
        elif p != 2 and alpha == 0 and kronecker(-1,p)==1:
            return -1
        else:
            return 1
    else:
        if p == 2 and r == 3 and alpha == 0:
            return 1
        else:
            return 0



def beta_rhopm_f(n,f):
    ret = 1
    for p,r in factor(n):
        ret *= beta_rhopm_f_local(p,r,vp(f,p))
    return ret



###########################################################

def c3(k):
    if k % 3 == 0:
        return -1/3
    elif k % 3 == 1:
        return 0
    else:
        return 1/3

def c4(k):
    if k % 4 == 0:
        return -1/4
    elif k % 4 == 1:
        return 0
    elif k % 4 == 2:
        return 1/4
    else:
        return 1/2

def c0(k,f):
    if k == 2 and f == 1:
        return 1
    else:
        return 0

def omega(n):
    return len(factor(n))

def beta(n):
    ret = 1
    for p,r in factor(n):
        if r == 1:
            ret *= -2
        elif r == 2:
            ret *= 1
        else:
            ret *= 0
    return ret 
    

###############################################



def get_dim_newspace(N,k,chi):
    f = chi.conductor()
    f_fact = factor(f)
    ret = 0
    ret += (k-1)/12 * psi(f) * beta_psi_f(N//f, f)
    ret += -1 * c3(k) * rho(f, chi, f_fact) * beta_rho_f(N//f, f)
    ret += -1 * c4(k) * rhopm(f, chi, f_fact) * beta_rhopm_f(N//f, f)
    ret += (-1/2) * 2^omega(f) * beta_sigma_f(N//f, f) 
    ret += c0(k,f) * moebius(N//f)
    return ret


################################################################


def test_small_dim(N_ub):
    for N in range(1,N_ub):
        for chi_idx, chi in enumerate(DirichletGroup(N)):
            dims = []
            for k in range(2,10):
                if chi(-1) == (-1)^k:
                    expdim = get_dim_newspace(N,k,chi)
                    assert expdim == dimension_new_cusp_forms(chi,k)
                    dims.append(expdim)
            sgn_ = '+' if chi(-1) == 1 else '-'
            print(f'({N},{chi_idx}){sgn_}: {dims}')



def test_k13(N_ub):
    for N in range(1,N_ub):
        for k in [13,25,37]:
            dims = {d:[] for d in divisors(N)}
            for chi_idx, chi in enumerate(DirichletGroup(N)):
                if chi(-1) == (-1)^k:
                    expdim = get_dim_newspace(N,k,chi)
                    assert expdim == dimension_new_cusp_forms(chi,k)
                    dims[chi.conductor()].append(expdim)
            print(f'({N},{k}) - ' + ''.join(f'{d}: {dims[d]},  ' for d in divisors(N) if len(dims[d])>0))




def compute_dims_smallNk(N_ub):
    for N in range(1, N_ub):
        for chi_idx, chi in enumerate(DirichletGroup(N)):
            ff = chi.conductor()
            if ff % 2 == 0 and N//ff % 4 == 2:
                continue
            if 1 <= N <= 50: kUB = 100
            elif 51 <= N <= 200: kUB = 20
            else: kUB = 8
            for k in range(2,kUB):
                if chi(-1) == (-1)^k:
                    expdim = get_dim_newspace(N,k,chi)
                    if expdim <= 1:
                        print(f'({N},{k},{chi.conrey_number()}): {expdim}')






def compute_dims(filenm):
    with open(filenm) as fl:
        ALL_N_kLB = eval(fl.read())
    for N, kLB in ALL_N_kLB:
        for chi_idx, chi in enumerate(DirichletGroup(N)):
            ff = chi.conductor()
            if ff % 2 == 0 and N//ff % 4 == 2:
                continue
            for k in range(2, kLB):
                if chi(-1) == (-1)^k:
                    expdim = get_dim_newspace(N,k,chi)
                    if expdim <= 1: 
                        print(f'({N},{k},{chi.conrey_number()}): {expdim}')




if ARG == 'test':
    test_small_dim(int(ARG_PARAM))
elif ARG == 'test_k13':
    test_k13(int(ARG_PARAM))
elif ARG == 'small_N_k':
    compute_dims_smallNk(int(ARG_PARAM))
elif ARG == 'N_kLB':
    compute_dims(ARG_PARAM)
else:
    assert False


