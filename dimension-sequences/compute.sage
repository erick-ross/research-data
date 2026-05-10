import sys
from sage.modular.dims import dimension_cusp_forms, dimension_new_cusp_forms


def __main__():
    ARG = sys.argv[1]
    IS_NS = {'FS': False, 'NS': True}[sys.argv[2]]
    if   ARG == 'search-level':
        search_level(is_NS=IS_NS)
    elif ARG == 'search-weight':
        search_weight(is_NS=IS_NS)
    elif ARG == 'search-sgnpatt-weight':
        search_sgnpatt_weight(is_NS=IS_NS)
    elif ARG == 'display-densities':
        display_densities(is_NS=IS_NS)
    elif ARG == 'display-density1-seqs':
        assert IS_NS == False, 'NS unnecessary here'
        display_density1_sgnpatt_seqs()
    elif ARG == 'verify-Tr-Wd-formula':
        N_UB = int(sys.argv[3])
        verify_Tr_Wd_formula(N_UB, is_NS=IS_NS)
    else:
        assert False, 'invalid parameter'



########## Search level-indexed sequences ##############################################


def B_k(k, N):
    ret = 0
    ret += (2*k-1)/12 * N
    ret -= (1/2) * 4.862 * N^(3/4)
    ret -= (5/6) * 4.862 * N^(1/4)
    return ret.n()

def B_k_pm(k, N):
    ret = 0
    ret += (2*k-1)/12 * (1/9.930) * N^(31/32)
    ret -= (1/2) * N^(1/2) 
    ret -= (5/6) * 4.862 * N^(1/4)
    ret -= 1
    return ret.n()



def search_level(is_NS=False):
    if is_NS: 
        dim_formula = dimension_new_cusp_forms
        dim_LB = B_k_pm
        N0_threshold = 4000
        k_vals = [7,6,5,4,3,2,1]
    else:
        dim_formula = dimension_cusp_forms
        dim_LB = B_k
        N0_threshold = 240000
        k_vals = [7,5,4,3,2,1]
    N_UB = 15000000  
    DIM_UB = 2000000

    for k in k_vals:
        print(f'RUNNING search_level(is_NS={is_NS}) for k={k}')
        dim_fnd = [0] * DIM_UB
        min_unfnd = 0
        for N in range(1, N_UB):
            dim = dim_formula(N, 2*k)
            dim_fnd[dim] = 1
            while dim_fnd[min_unfnd] == 1:
                min_unfnd += 1
            if N >= N0_threshold and dim_LB(k, N) > min_unfnd:
                # We have already computed all possible dims <= dim_LB(k, N)
                # And min_unfnd was not one of them
                print(f'FOUND: {min_unfnd} is not a dim for k={k}')
                print(f'N0 = {N}, Bk({N}) = {dim_LB(k, N)}')
                break
    print(f'DONE!!!')
    
   
####### Search weight-indexed sequences #########################################################

def Lemma_4_1(a,b,seq):
    fnd_values = [0]*(2*a)
    for k in range(1, 24*b+2):
        if 0 <= seq[k] < 2*a:
            fnd_values[ seq[k] ] = 1
    return all(fnd_values)


def psi(N):
    ret = 1
    for p,r in factor(N):
        ret *= (p+1)*p^(r-1)
    return ret

def psi_pm(N):
    ret = 1
    for p,r in factor(N):
        if r == 1:   ret *= p-1
        elif r == 2: ret *= p^2-p-1
        else:        ret *= (p^2-1)*(p-1)*p^(r-3)
    return ret

def omega(N):
    return len(factor(N))

def search_weight(is_NS=False):
    print(f'RUNNING search_weight(is_NS={is_NS})')
    if is_NS:
        dim_formula = dimension_new_cusp_forms
        N_vals = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 16, 18, 22, 25, 28, 60]
        psi_func = psi_pm
    else:
        dim_formula = dimension_cusp_forms
        N_vals = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 16, 18, 25]
        psi_func = psi

    good_seqs = []
    for N in N_vals:
        a = 2 * psi_func(N)
        b = 1
        seq = [None] + [dim_formula(N,2*k) for k in range(1,24*b+2)]
        all_nums = Lemma_4_1(a, b, seq)
        print(f'N={N} : all_nums={all_nums}')
        print(seq)
        if all_nums:
            good_seqs.append(N)
    print('good sequences:', good_seqs)




######### Implementation of Atkin-Lehner trace formula ####################################

# In this section, we implement an explicit Tr Wd formula from Kimball Martin.
# The periodicity of Tr Wd (for d != 1, Tr Wd is 12-periodic for k >= 2) is clear from this formula.
# This is Equation 1.6 of "Refined dimensions of cusp forms, and equidistribution and bias of signs"
# We use this explicit formula because computing Tr Wd natively in Sage is too slow.
# To match Martin's notation, we use the variables:  
# kk (which is our 2*k), M (which is != 1, and is our d), and M_pm (which is our N/d)


def get_a(M, M_pm):
    if M % 8 in [1,2,5,6]:
        return 1
    elif M % 8 == 3:
        if M_pm % 2 == 1: return 4
        else:             return 6
    elif M % 8 == 7:
        if M_pm % 2 == 1: return 2
        else:             return 4

# There are more efficient ways to do this.
# But this works fine for our purposes
def get_r(D, n):
    ret = 0
    for r in range(2*n):
        if (r^2 - D) % (4*n) == 0:
            ret += 1
    return ret

def odd_part(M_pm):
    while M_pm % 2 == 0:
        M_pm //= 2
    return M_pm

def p_kk_sqrt2(kk):
    if   kk % 8 == 0: return -1
    elif kk % 8 == 2: return 1
    elif kk % 8 == 4: return 1
    elif kk % 8 == 6: return -1

def p_kk_sqrt3(kk):
    if   kk % 12 == 0:  return -1
    elif kk % 12 == 2:  return 1
    elif kk % 12 == 4:  return 2
    elif kk % 12 == 6:  return 1
    elif kk % 12 == 8:  return -1
    elif kk % 12 == 10: return -2
        
# Data from https://oeis.org/A014600 (run 'sage parse_oeis_data.sage < oeis_data.txt')
H_PM = {
    -4: 1/2, -8: 1, -3: 1/3, -20: 2, -24: 2, -7: 1, -40: 2, -11: 1, -52: 2, -56: 4, 
    -15: 2, -84: 4, -88: 2, -104: 6, -120: 4, -132: 4, -68: 4, -136: 4, -35: 2, -19: 1, 
    -152: 6, -39: 4, -168: 4, -23: 3, -184: 4, -264: 8, -280: 4, -312: 4, -51: 2, 
    -408: 4, -420: 8, -55: 4, -440: 12, -228: 4, -456: 8, -260: 8, -520: 4, -276: 8, 
    -552: 8, -840: 8, -660: 8, -1320: 8, -195: 4, -1560: 16, -116: 6, -31: 3, -148: 2
}

def delta_M(M):
    if (-M) % 4 == 1: return -M
    else:             return -4*M
    

def Tr_WM_FS_formula(N, kk, M):
    assert is_squarefree(N)
    assert M != 1
    M_pm = N // M
    delta = delta_M(M)
    M_pm_odd = M_pm
    while M_pm_odd % 2 == 0: M_pm_odd //= 2
    ret = 0
    ret += (1/2) * (-1)^(kk//2) * get_a(M, M_pm) * H_PM[delta] * get_r(delta, M_pm_odd)
    ret -= (1/2) * kronecker_delta(M,2) * p_kk_sqrt2(kk) * get_r(-4, M_pm)
    ret -= (1/3) * kronecker_delta(M,3) * p_kk_sqrt3(kk) * get_r(-3, M_pm)
    ret += kronecker_delta(kk,2)
    return ZZ(ret)

def Tr_WM_NS_formula(N, kk, M):
    assert is_squarefree(N)
    assert M != 1
    M_pm = N // M
    ret = 0
    for tt in divisors(M_pm):
        ret += (-2)^omega(M_pm//tt) * Tr_WM_FS_formula(tt*M, kk, M)
    return ret


# The end-user function we will use (with variables k and d)
def get_Tr_Wd(N, k, d, is_NS=False):
    if is_NS:
        if d == 1: return dimension_new_cusp_forms(N, 2*k)
        else:      return Tr_WM_NS_formula(N, 2*k, d)
    else:     
        if d == 1: return dimension_cusp_forms(N, 2*k)
        else:      return Tr_WM_FS_formula(N, 2*k, d) 


# Verify the explicit formula from Kimball Martin against native Sage computations
# We only need to check 1 <= k < 14 because Tr Wd is 12-periodic for k >= 2
def verify_Tr_Wd_formula(N_UB, is_NS=False):
    print(f'RUNNING verify_Tr_Wd_formula({N_UB}, is_NS={is_NS})')
    k_UB = 14
    for N in range(1, N_UB):
        if not is_squarefree(N):
            continue
        for k in range(1, k_UB):
            S = ModularSymbols(N, 2*k, sign=1).cuspidal_submodule()
            if is_NS: S = S.new_submodule()
            print_lst = []
            for d in divisors(N):
                mat = d^(-(k-1)) * S.atkin_lehner_operator(d).matrix()
                sage_trace = ZZ(mat.trace())
                print_lst.append(sage_trace)
                assert sage_trace == get_Tr_Wd(N, k, d, is_NS=is_NS)
            print(f'N={N}, k={k}: {print_lst}')
    print('DONE!!')





##### Search sign-pattern weight indexed sequences ##########################################

def get_sgnpatts(N):
    assert is_squarefree(N)
    pms = [p for p,_ in factor(N)]
    t = len(pms)
    ret = []
    for mask in range(2^t):
        sgnpatt = {}
        for i in range(t):
            sgnpatt[ pms[i] ] = [1, -1][ (mask>>i)&1 ]
        ret.append(sgnpatt)
    return ret

def eval_sgnpatt(sgnpatt, t):
    assert is_squarefree(t)
    ret = 1
    for p,r in factor(t):
        ret *= sgnpatt[p]
    return ret

def sgnpatt_str(sgnpatt):
    pms = sorted(sgnpatt.keys())
    return ''.join({1: '+', -1: '-'}[sgnpatt[p]]  for p in pms)
    

def construct_sgnpatt_dim_table(N, k_UB, is_NS=False):
    if is_NS: dim_formula = dimension_new_cusp_forms
    else:     dim_formula = dimension_cusp_forms

    table = {}
    for sgnpatt in get_sgnpatts(N):
        dims = [None]
        for k in range(1, k_UB):
            dim = 0
            for d in divisors(N):
                dim += eval_sgnpatt(sgnpatt, d) * get_Tr_Wd(N, k, d, is_NS=is_NS)
            assert dim % 2^omega(N) == 0
            dims.append( dim // 2^omega(N) )
        table[sgnpatt_str(sgnpatt)] = dims
    # gut check that the total dimensions are correct
    for k in range(1, k_UB):
        total_dim =  sum(table[sp_str][k] for sp_str in table)
        assert dim_formula(N, 2*k) == total_dim
    return table




def search_sgnpatt_weight(is_NS=False):
    print(f'RUNNING search_sgnpatt_weight(is_NS={is_NS})')
    if is_NS:
        N_vals = [1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 21, 22, 26, 30, 33, 34, 35, 38, 
                    39, 42, 46, 66, 70, 78, 102, 105, 110, 114, 130, 138, 210, 330, 390]
        psi_func = psi_pm
    else:
        N_vals = [1, 2, 3, 5, 6, 7, 10, 11, 14, 15]
        psi_func = psi

    good_seqs = []
    for N in N_vals:
        a = 2 * psi_func(N)
        b = 2^omega(N)
        sgnpatt_dim_table = construct_sgnpatt_dim_table(N, 24*b+2, is_NS=is_NS)
        for sp_str in sgnpatt_dim_table:
            all_nums = Lemma_4_1(a, b, sgnpatt_dim_table[sp_str])
            print(f'N = {N}, sgnpatt = {sp_str} : all_nums={all_nums}')
            print(sgnpatt_dim_table[sp_str])
            if all_nums:
                good_seqs.append((N,sp_str))
    print('good sequences:', good_seqs)



#### Display the density properties noted in the paper ##########################################


def display_densities(is_NS=False):
    print(f'RUNNING display_densities(is_NS={is_NS})')
    if is_NS: psi_func = psi_pm
    else:     psi_func = psi
    dens_ge1 = []
    dens_eq1 = []
    for N in range(1, 1000000):
        if is_squarefree(N):
            dens = 6 * 2^omega(N) / psi_func(N)
            if dens >= 1:
                print(N, dens)
                dens_ge1.append(N)
                if dens == 1:
                    dens_eq1.append(N)
    print(dens_ge1)
    print(dens_eq1)




def display_density1_sgnpatt_seqs():
    print(f'RUNNING display_density1_sgnpatt_seqs()')
    N_vals = [11,14,15]
    for N in N_vals:
        a = 2 * psi(N)
        b = 2^omega(N)
        sgnpatt_dim_table = construct_sgnpatt_dim_table(N, 14, is_NS=False)
        for sp_str in sgnpatt_dim_table:
            print(f'N = {N}, sgnpatt = {sp_str} :')
            print(sgnpatt_dim_table[sp_str])

##########################################################################            
__main__()


