import sys 

ARG = sys.argv[1]
if ARG == 'test':
    M_VALUES = list(range(1,17))
    K_MIN, K_MAX, N_MAX = 2, 16, 20
elif ARG == 'check_decr':
    M_VALUES = [1,2,4]
    K_MIN, K_MAX, N_MAX = 2, 2000, 3392663
elif ARG == 'check_nonrep_A':
    M_VALUES = [1,4,16]
    K_MIN, K_MAX, N_MAX = 2, 18000, 100
elif ARG == 'check_nonrep_B':
    M_VALUES = [1,4,16]
    K_MIN, K_MAX, N_MAX = 2, 200, 332427
elif ARG == 'check_T2_T3':
    M_VALUES = [1,2,3,4,9]
    K_MIN, K_MAX, N_MAX = 2, 200, 1
elif ARG == 'check_conj':
    M_VALUES = list(range(1,11)) + [16, 25, 36, 49, 64, 81, 100]
    K_MIN, K_MAX, N_MAX = 2, 200, 100
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


_6_h_w = {
         -3: 2, -4: 3, -7: 6, -8: 6, -11: 6, -12: 6, -15: 12, -16: 6, -19: 6, -20: 12,
         -23: 18, -24: 12, -27: 6, -28: 6, -31: 18, -32: 12, -35: 12, -36: 12, -39: 24, 
         -40: 12, -43: 6, -44: 18, -47: 30, -48: 12, -51: 12, -52: 12, -55: 24, -56: 24, 
         -59: 18, -60: 12, -63: 24, -64: 12, -67: 6, -68: 24, -71: 42, -72: 12, -75: 12, 
         -76: 18, -79: 30, -80: 24, -83: 18, -84: 24, -87: 36, -88: 12, -91: 12, -92: 18, 
         -95: 48, -96: 24, -99: 12, -100: 12, -103: 30, -104: 36, -107: 18, -108: 18, 
         -111: 48, -112: 12, -115: 12, -116: 36, -119: 60, -120: 24, -123: 12, -124: 18, 
         -127: 30, -128: 24, -131: 30, -132: 24, -135: 36, -136: 24, -139: 18, -140: 36, 
         -143: 60, -144: 24, -147: 12, -148: 12, -151: 42, -152: 36, -155: 24, -156: 24, 
         -159: 60, -160: 24, -163: 6, -164: 48, -167: 66, -168: 24, -171: 24, -172: 18, 
         -175: 36, -176: 36, -179: 30, -180: 24, -183: 48, -184: 24, -187: 12, -188: 30, 
         -191: 78, -192: 24, -195: 24, -196: 24, -199: 54, -200: 36, -203: 24, -204: 36, 
         -207: 36, -208: 24, -211: 18, -212: 36, -215: 84, -216: 36, -219: 24, -220: 24, 
         -223: 42, -224: 48, -227: 30, -228: 24, -231: 72, -232: 12, -235: 12, -236: 54, 
         -239: 90, -240: 24, -243: 18, -244: 36, -247: 36, -248: 48, -251: 42, -252: 24, 
         -255: 72, -256: 24, -259: 24, -260: 48, -263: 78, -264: 48, -267: 12, -268: 18, 
         -271: 66, -272: 48, -275: 24, -276: 48, -279: 72, -280: 24, -283: 18, -284: 42, 
         -287: 84, -288: 24, -291: 24, -292: 24, -295: 48, -296: 60, -299: 48, -300: 36, 
         -303: 60, -304: 36, -307: 18, -308: 48, -311: 114, -312: 24, -315: 24, -316: 30, 
         -319: 60, -320: 48, -323: 24, -324: 36, -327: 72, -328: 24, -331: 18, -332: 54, 
         -335: 108, -336: 48, -339: 36, -340: 24, -343: 42, -344: 60, -347: 30, -348: 36, 
         -351: 72, -352: 24, -355: 24, -356: 72, -359: 114, -360: 48, -363: 24, -364: 36, 
         -367: 54, -368: 36, -371: 48, -372: 24, -375: 60, -376: 48, -379: 18, -380: 48, 
         -383: 102, -384: 48, -387: 24, -388: 24, -391: 84, -392: 48, -395: 48, -396: 36, 
         -399: 96, -400: 24
        }




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



print('computing MU_SUM values')
MU_SUM = {}
for m in M_VALUES:
    for t,_ in get_t_2wt(m):
        for n in get_n(m, t):
            # mu_sum is a multiplicative function of N; see Cohen+Stromberg Remark 12.4.12
            mu_sm = [0]*(N_MAX+1)
            for N in range(1, N_MAX+1):
                N_fact = factor(N)
                if len(N_fact) == 1:
                    mu_sm[N] = mu_sum(N,t,n,m)
                else:
                    mu_sm[N] = product(mu_sm[pm^exp] for (pm,exp) in N_fact)
            MU_SUM[(t,n,m)] = mu_sm



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


def get_dim(N, k):
    return get_trace(1, N, k)

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



########################################################################################
# The following eight functions make up the only section of the code specific to the paper


def get_T2_E(N):
    N_fact = factor(N)
    _2_omega_N = 2^len(N_fact)
    sigma0_N = product([r+1 for (p,r) in N_fact])
    psi_N = product([(p+1)*p^(r-1) for (p,r) in N_fact])
    theta1 = _2_omega_N * sqrt(N) / psi_N
    theta2 = sigma0_N / psi_N
    theta3 = _2_omega_N^2 / psi_N
    theta4 = _2_omega_N / psi_N
    theta5 = 1 / psi_N
    return (7/2)*theta1 + theta2 + (129/2+8*sqrt(2))*theta3 + (24*sqrt(2)+323/6)*theta4 + 9*theta5


def get_T2_k_N(N):
    E_N = get_T2_E(N).n()
    k_N = 1
    while (6*k_N+5)/24 <= E_N:
        k_N += 1
    return k_N



def check_decreasing(N):
    k_N = get_T2_k_N(N)
    if N < 1000 or N % 10000 == 1:
        print(f'Checking N={N} ; \t\t k_N={k_N}' )

    # Now verify that a2(T2(N,2*k)) is strictly decreasing 
    # for 1 <= k <= k_N  s.t. s(N,2k) >= 2
    if k_N > 1:
        a2_values = []
        for k in range(1, k_N+1):
            if get_dim(N, 2*k) >= 2:
                a2_values.append(get_a2_coeff(2, N, 2*k))
        for i in range(len(a2_values)-1):
            assert a2_values[i] > a2_values[i+1]
        


def get_T4_E(N):
    N_fact = factor(N)
    _2_omega_N = 2^len(N_fact)
    psi_N = product([(p+1)*p^(r-1) for (p,r) in N_fact])
    theta1 = _2_omega_N * sqrt(N) / psi_N
    theta4 = _2_omega_N / psi_N
    theta5 = 1 / psi_N
    val = (17/2) * theta4 + (1/4) * theta1
    return (17/12)*val + (1/3)*val^2 + (234577/576)*theta5


def get_T4_k_N(N):
    E_N = get_T4_E(N).n()
    k_N = 1
    while (10*k_N+5)/192 <= E_N:
        k_N += 1
    return k_N



def check_nonrep(N):
    k_N = get_T4_k_N(N)
    if N < 1000 or N % 10000 == 1:
        print(f'Checking N={N} ; \t\t k_N={k_N}' )

    # Now verify that each a2(T4(N,2k)) < a2(T4(N,2k_N)) 
    # for 1 <= k < k_N  s.t. s(N,2k) >= 2, 
    # AND that these values are non-repeating. 
    if k_N > 1:
        a2_values = []
        for k in range(1, k_N):
            if get_dim(N, 2*k) >= 2:
                a2_values.append(get_a2_coeff(4, N, 2*k))
        a2_k_N = get_a2_coeff(4, N, 2*k_N)
        assert max(a2_values) < a2_k_N
        assert len(a2_values) == len(set(a2_values))
        



def check_T2_T3():
    for k in [12] + list(range(14, 82)):
        a2_T2 = get_a2_coeff(2, 1, 2*k)
        a2_T3 = get_a2_coeff(3, 1, 2*k)
        assert a2_T2 != a2_T3



def check_conj():
    for m in range(2,11):
        for N in range(1, N_MAX+1):
            if gcd(m,N) > 1:
                continue
            print(f'checking conjecture at m={m}, N={N}')
            a2_values = []
            for k in range(1, K_MAX//2+1):
                if get_dim(N, 2*k) >= 2:
                    a2_values.append(get_a2_coeff(m, N, 2*k))
            
            assert len(a2_values) == len(set(a2_values))




##############################################################################




if ARG == 'test':
    test_all_traces()
    for m in range(1,5):
        test_a2_coeff(m)
elif ARG == 'check_decr':
    for N in range(1, N_MAX+1, 2):
        check_decreasing(N)
    print(f'FINISHED: For each odd N <= {N_MAX}, a2(T2(N,2k)) is strictly decreasing in k')
elif ARG == 'check_nonrep_A':
    for N in range(1, 100, 2):
        check_nonrep(N)
    print(f'FINISHED: For each odd 1 <= N <= 99, a2(T4(N,2k)) is non-repeating in k')
elif ARG == 'check_nonrep_B':
    for N in range(101, N_MAX+1, 2):
        check_nonrep(N)
    print(f'FINISHED: For each odd 101 <= N <= {N_MAX}, a2(T4(N,2k)) is non-repeating in k')
elif ARG == 'check_T2_T3':
    check_T2_T3()
    print('FINISHED: a2(T2(1,2k)) != a2(T3(1,2k)) for 12 <= k <= 81, k != 13')
elif ARG == 'check_conj':
    check_conj()


