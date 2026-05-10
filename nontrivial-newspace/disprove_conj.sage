import sys



#######################################3


def beta_psi_f_local(p,r):
    assert r >= 1
    if r == 1:
        return p-1
    elif r==2:
        return p^2-p-1
    else:
        return (p^3-p^2-p+1) * p^(r-3)


def beta_psi_f(n):
    ret = 1
    for p,r in factor(n):
        ret *= beta_psi_f_local(p,r)
    return ret


################################################

def beta_sigma_f_local(p,r):
    assert r >= 1
    if r % 2 == 1:
        return 0
    elif r == 2:
        return p-2
    else:
        return (p^2-2*p+1) * p^(r//2-2)
    


def beta_sigma_f(n):
    ret = 1
    for p,r in factor(n):
        ret *= beta_sigma_f_local(p,r)
    return ret

#######################################################


def beta_rho_f_local(p,r):
    assert r >= 1
    if r == 1:
        if p == 3:
            return -1
        elif p != 2 and p != 3 and kronecker(-3,p)==1:
            return 0
        else:
            return -2
    elif r == 2:
        if p == 3:
            return -1
        elif p != 2 and p != 3 and kronecker(-3,p)==1:
            return -1
        else:
            return 1
    else:
        if p == 3 and r == 3:
            return 1
        else:
            return 0



def beta_rho_f(n):
    ret = 1
    for p,r in factor(n):
        ret *= beta_rho_f_local(p,r)
    return ret



##################################################3



def beta_rhopm_f_local(p,r):
    assert r >= 1
    if r == 1:
        if p == 2:
            return -1
        elif p != 2 and kronecker(-1,p)==1:
            return 0
        else:
            return -2
    elif r == 2:
        if p == 2:
            return -1
        elif p != 2 and kronecker(-1,p)==1:
            return -1
        else:
            return 1
    else:
        if p == 2 and r == 3:
            return 1
        else:
            return 0



def beta_rhopm_f(n):
    ret = 1
    for p,r in factor(n):
        ret *= beta_rhopm_f_local(p,r)
    return ret



###########################################################





def get_dim_newspace(N):
    ret = 0
    ret += 1/12 * beta_psi_f(N)
    ret += -1/3 * beta_rho_f(N)
    ret += -1/4 * beta_rhopm_f(N)
    ret += -1/2 * beta_sigma_f(N) 
    ret += moebius(N)
    return ret


################################################################


def test_small_dim(N_ub):
    for N in range(1,N_ub):
        expdim = get_dim_newspace(N)
        assert expdim == dimension_new_cusp_forms(N,2)
        print(f'{N}: {expdim}')





def disprove_conj():
    UB = 58260767
    dim_UB = 67846

    cnt = [0] * 10000000
    for N in range(1, UB):
        dim = get_dim_newspace(N)
        cnt[dim] += 1

    
    for i in range(dim_UB+1):
        if cnt[i] == 0:
            print('NOT FOUND:', i)
    
    
    


# test_small_dim(50000)
disprove_conj()









