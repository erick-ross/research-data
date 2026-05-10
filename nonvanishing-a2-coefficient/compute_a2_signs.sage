import time, sys 

ARG = sys.argv[1]
if ARG == 'test':
    M_VALUES = list(range(1,17))
    K_MIN, K_MAX, N_MAX = 2, 16, 20
elif ARG[0] == '2':
    M_VALUES = [1,2,4]
    BOUNDS = {
        '2_2M':   (2, 4, 2700000-1),
        '2_150K': (6, 14, 150000-1),
        '2_8K':   (16, 48, 8800-1),
        '2_571':  (50, 146, 571-1),
        '2_43':   (148, 560, 43-1),
        '2_small': (2, 40, 5000)
    }
    K_MIN, K_MAX, N_MAX = BOUNDS[ARG]
elif ARG[0] == '3':
    M_VALUES = [1,3,9]
    BOUNDS = {
        '3_63M':  (2, 2, 63000000-1),
        '3_2M':   (4, 8, 2700000-1),
        '3_150K': (10, 32, 150000-1),
        '3_8K':   (34, 114, 8800-1),
        '3_571':  (116, 344, 571-1),
        '3_43':   (346, 1288, 43-1),
        '3_small': (2, 40, 5000)
    }
    K_MIN, K_MAX, N_MAX = BOUNDS[ARG]
elif ARG[0] == '4':
    M_VALUES = [1,4,16]
    BOUNDS = {
        '4_2M':   (2, 2, 2700000-1),
        '4_150K': (4, 12, 150000-1),
        '4_8K':   (14, 122, 8800-1),
        '4_571':  (124, 1496, 571-1),
        '4_43':   (1498, 62940, 43-1),
        '4_small': (2, 200, 2000)
    }
    K_MIN, K_MAX, N_MAX = BOUNDS[ARG]
else:
    assert False, 'invalid argument: ' + ARG 




##############################################

def psi(N):
    ret = 1
    for pm, exp in factor(N):
        ret *= (pm+1) * pm^(exp-1)
    return ret

print('computing PSI values')
PSI = [0]*(N_MAX+1)
for N in range(1,N_MAX+1):
    PSI[N] = psi(N)

def _12_A1(m, N, k):
    sqrtm = round(sqrt(m))
    if sqrtm^2 != m or gcd(sqrtm, N) != 1:
        return 0
    else:
        return (k-1) * PSI[N] * m^(k//2 - 1)
    

####################################################

def get_t_2wt(m):
    ret = []
    for t in range(4*m+1):
        if t^2 >= 4*m: 
            break
        _2wt = 1 if (t == 0) else 2
        ret.append((t, _2wt))
    return ret


def get_n(m,t):
    ret = []
    for n in range(1, 4*m - t^2 + 1):
        if (t^2 - 4*m) % (n^2) == 0 and ((t^2 - 4*m)//(n^2))%4 in [0,1]:
            ret.append(n)
    return ret




U = {}
for m in M_VALUES:
    for t,_ in get_t_2wt(m):
        # create U_table
        U_ = [0]*(K_MAX+1)
        U_[0] = 0
        U_[1] = 1
        for i in range(2,K_MAX+1):
            U_[i] = t * U_[i-1] - m * U_[i-2]
        U[(t,m)] = U_

_6_h_w = {-3: 2, -4: 3, -7: 6, -8: 6, -11: 6, -12: 6, -15: 12, -16: 6, -19: 6, -20: 12, 
        -23: 18, -24: 12, -27: 6, -28: 6, -31: 18, -32: 12, -35: 12, -36: 12, -39: 24,
        -40: 12, -43: 6, -44: 18, -47: 30, -48: 12, -51: 12, -52: 12, -55: 24, -56: 24,
        -59: 18, -60: 12, -63: 24, -64: 12, -67: 6, -68: 24, -71: 42, -72: 12, -75: 12, 
        -76: 18, -79: 30, -80: 24, -83: 18, -84: 24, -87: 36, -88: 12, -91: 12, -92: 18, 
        -95: 48, -96: 24, -99: 12, -100: 12, -103: 30, -104: 36, -107: 18, -108: 18, -111: 48,
        -112: 12, -115: 12, -116: 36, -119: 60, -120: 24, -123: 12, -124: 18, -127: 30, 
        -128: 24, -131: 30, -132: 24, -135: 36, -136: 24, -139: 18, -140: 36, -143: 60, -144: 24}


def mu_sum(N,t,n,m):
    R.<x>=PolynomialRing(Integers(N))
    Nn = gcd(N,n)
    ret = 0
    poly = x^2-t*x+m 
    roots_modN = poly.roots(multiplicities=False)
    for soln_ in roots_modN:
        soln = int(soln_)
        if gcd(soln, N) != 1:
            continue
        # check if it solves the eqn for some lifting of soln
        soln_lifts = False
        for i in range(Nn):
            soln_lftd = soln + i*N
            value = soln_lftd^2 - t*soln_lftd + m
            if value % (N*Nn) == 0:
                soln_lifts = True
                break
        if soln_lifts:
            ret += 1
    return ret




MU_SUM = {}
for m in M_VALUES:
    for t,_ in get_t_2wt(m):
        for n in get_n(m, t):
            t0 = time.time()
            print('computing MU_SUM', m,t,n)
            # mu_sum is a multiplicative function of N; see Cohen+Stromberg Remark 12.4.12
            mu_sm = [0]*(N_MAX+1)
            for N in range(1, N_MAX+1):
                N_fact = factor(N)
                if len(N_fact) == 1:
                    mu_sm[N] = mu_sum(N,t,n,m)
                else:
                    mu_sm[N] = product(mu_sm[pm^exp] for (pm,exp) in N_fact)
            MU_SUM[(t,n,m)] = mu_sm
            print(time.time() - t0)



def mu(N,t,n,m):
    Nn = gcd(N,n)
    return (PSI[N] // PSI[N//Nn]) * MU_SUM[(t,n,m)][N]


def _12_A2(m,N,k):
    ret = 0 
    for t, _2wt in get_t_2wt(m):
        for n in get_n(m,t):
            ret += _2wt * U[(t,m)][k-1] * _6_h_w[(t^2 - 4*m)//(n^2)] * mu(N,t,n,m)
    return ret



######################################################################


# d <= sqrt(m) with weight 2 if d != sqrt(m)
def get_d_2wt(m):
    ret = []
    for d in divisors(m):
        if d^2 < m:
            ret.append((d,2))
        elif d^2 == m:
            ret.append((d,1))
    return ret



def _12_A3(m,N,k):
    ret = 0
    for d,_2_d_wt in get_d_2wt(m):
        for tau in divisors(N):
            g1 = gcd(tau, N//tau)
            g2 = gcd(N, m//d - d)
            if g2 % g1 != 0: 
                continue
            y = CRT([d,m//d], [tau,N//tau])
            if gcd(y, N) > 1:
                continue
            ret += _2_d_wt * d^(k-1) * euler_phi(gcd(tau, N//tau)) 
    return 6*ret
    

def _12_A4(m,N,k):
    if k > 2:
        return 0
    ret = 0
    for c in divisors(m):
        if gcd(N, m//c) == 1:
            ret += c
    return 12 * ret


######################################################################

def get_trace(m,N,k):
    A1 = _12_A1(m,N,k) 
    A2 = _12_A2(m,N,k) 
    A3 = _12_A3(m,N,k) 
    A4 = _12_A4(m,N,k)
    ret = A1 - A2 - A3 + A4
    # print(N,m,k, '|', A1, A2, A3, A4, '|',  ret//12)
    assert ret % 12 == 0
    return ret // 12


def get_a2_coeff(m, N, k):
    ret = get_trace(m, N, k)^2
    for d in divisors(m):
        if gcd(d, N) == 1:
            ret -= d^(k-1) * get_trace((m//d)^2, N, k)
    assert ret % 2 == 0
    return ret // 2



def test_all_traces():
    for N in range(1, N_MAX+1):
        for m in M_VALUES:
            for k in range(K_MIN, K_MAX+1, 2):
                tr = get_trace(m,N,k)
                print(N,m,k, 'Tr', tr)
                S = ModularForms(Gamma0(N),k).cuspidal_subspace()
                MS = ModularSymbols(Gamma0(N),k,sign=1).cuspidal_subspace()
                assert tr == S.hecke_matrix(m).trace()
                assert tr == MS.hecke_operator(m).trace()


def test_a2_coeff(m):
    for N in range(1, N_MAX+1):
        for k in range(K_MIN, K_MAX+1, 2):
            a2_coeff = get_a2_coeff(m,N,k)
            print(N,m,k, 'a2', a2_coeff)
            MS = ModularSymbols(Gamma0(N),k,sign=1).cuspidal_subspace()
            cply = [0,0,0] + MS.hecke_operator(m).charpoly().list()
            assert a2_coeff == cply[-3]



def find_all_nonneg_a2(m):
    N100 = max(2, N_MAX//100) 
    t0 = time.time()
    for N in range(1, N_MAX+1):
        if (N-1) % N100 == 0:
            print('computing N =', N,  time.time() - t0)
            t0 = time.time()
        if gcd(N,m) != 1:
            continue
        for k in range(K_MIN, K_MAX+1, 2):
            a2_coeff = get_a2_coeff(m,N,k)
            if a2_coeff >= 0:
                print('Non-Neg a2:  (m,N,k):', m, N, k, '\t dim', get_trace(1,N,k), ' \ta2:', a2_coeff)



def find_all_nonpos_a2(m):
    N100 = max(2, N_MAX//100) 
    t0 = time.time()
    for N in range(1, N_MAX+1):
        if (N-1) % N100 == 0:
            print('computing N =', N,  time.time() - t0)
            t0 = time.time()
        if gcd(N,m) != 1:
            continue
        for k in range(K_MIN, K_MAX+1, 2):
            a2_coeff = get_a2_coeff(m,N,k)
            if a2_coeff <= 0:
                print('Non-Pos a2:  (m,N,k):', m, N, k, '\t dim', get_trace(1,N,k), ' \ta2:', a2_coeff)


##############################################################################

if ARG == 'test':
    test_all_traces()
    for m in range(1,5):
        test_a2_coeff(m)
elif ARG[0] == '2':
    find_all_nonneg_a2(2)
elif ARG[0] == '3':
    find_all_nonneg_a2(3)
elif ARG[0] == '4':
    find_all_nonpos_a2(4)


