import sys
##########################################################
# Implementation of Skoruppa-Zagier-Assaf trace formula
##########################################################


def main():
    ARG = sys.argv[1]
    if ARG == 'search-weight':
        m = ZZ(sys.argv[2])
        N = ZZ(sys.argv[3])
        sgnpatt = signpattern(sys.argv[4])
        if sys.argv[5] == 'all':
            k_ub = get_k_ub(N, m)
        else:
            k_ub = ZZ(sys.argv[5])
        search_weight_c2_Tm_sgnpatt(m, N, sgnpatt, k_ub)
    elif ARG == 'test_trace':
        ARG_N_ub = ZZ(sys.argv[2])
        ARG_is_NS = {'NS': True, 'FS': False}[sys.argv[3]]
        test_Tr_Tm_WQ(ARG_N_ub, is_NS=ARG_is_NS)
        test_Tr_Tm_WQ(ARG_N_ub, is_NS=ARG_is_NS)
    elif ARG == 'check_conj_c2':
        ARG_k_ub = ZZ(sys.argv[2])
        ARG_N_ub = ZZ(sys.argv[3])
        ARG_m_ub = ZZ(sys.argv[4])
        ARG_is_NS = {'NS': True, 'FS': False}[sys.argv[5]]
        check_conj_c2(ARG_k_ub, ARG_N_ub, ARG_m_ub, is_NS=ARG_is_NS)
    else:
        assert False, 'invalid parameter'




# This returns the values of the Lucas Sequence of the first kind. 
# It is only implemented for n odd, using the identity  U(2n+1) = U*(n)
# We use this identity because U* is an integer sequence
# Note the Q here is unrelated to the Q elsewhere
def get_U(P, Q, n):
    assert n % 2 == 1
    return Ustar_seq.seq(ZZ(P^2), Q).get(n//2)
# Here,  U*(n) = (P^2-2Q) U*(n-1) - Q^2 U(n-2)  and   U*(0)=1,  U*(1)=P^2-Q 
# This is most efficient when the values for n are accessed sequentially.
class Ustar_seq:
    SEQS = {}
    def seq(P2, Q):
        if (P2, Q) not in Ustar_seq.SEQS:
            Ustar_seq.SEQS[ (P2, Q)  ] = Ustar_seq(P2, Q)
        return Ustar_seq.SEQS[ (P2, Q) ]

    def __init__(self, P2, Q):
        self.P2 = P2
        self.Q = Q
        self.d = [1,P2-Q]
        self.cur = 1

    def get(self,n):
        if n < self.cur:
            self.d = [1,self.P2-self.Q]
            self.cur = 1
        # build up to n
        for i in range(self.cur+1, n+1):
            self.d[i%2] = (self.P2-2*self.Q) * self.d[(i-1)%2]  -  self.Q^2 * self.d[(i-2)%2]
            self.cur += 1
        return self.d[n%2]


    



def HK_classnum(D):
    return QQ(pari(abs(D)).qfbhclassno())


def sqfr_part(n):
    return product(p for (p,r) in factor(n) if r&1)



@cached_function
def H_N(N, Delta):
    a2b = gcd(N, abs(Delta))
    b = sqfr_part(a2b)
    a2 = a2b // b
    if abs(Delta) % (a2*b^2) != 0:
        return 0
    return a2b * kronecker(Delta//(a2*b^2), N//a2b) * HK_classnum(abs(Delta)//(a2*b^2))



@cached_function
def get_s(m,Q,Q_pm):
    ret = []
    for s in range(0, 4*m*Q_pm+1, Q_pm):
        if s^2 > 4*m*Q_pm:
            break
        vl = gcd((s//Q_pm)^2, Q//Q_pm)
        if not is_squarefree(vl):
            continue
        ret.append(s)
        if s != 0: 
            ret.append(-s)
    return ret




def get_B(N):
    ret2 = N // sqfr_part(N)
    return ZZ(sqrt(ret2))
        


def s_k_N(k,N,m,Q):
    term1 = 0
    for Q_pm in divisors(Q):
        for s in get_s(m,Q,Q_pm):
            U_k1 = get_U(s/sqrt(Q_pm), m, k-1) 
            term1 += U_k1 * H_N(N//Q, s^2-4*m*Q_pm)
        
    term2 = 0
    for m_pm in divisors(m):
        mn = min(m_pm, m//m_pm)
        B1 = gcd(get_B(Q), m_pm + m//m_pm)
        B2 = gcd(get_B(N//Q), m_pm - m//m_pm)
        term2 += mn^(k-1) * B1 * B2
    
    term3 = 0
    if k == 2 and is_square(N//Q):
        term3 = sigma(Q, 0) * sigma(m, 1)

    ret = (-1/2) * term1 + (-1/2) * term2 + term3
    return ZZ(ret)



def get_alpha(N):
    ret = 1
    for p,r in factor(N):
        if   r == 1: ret *= -1
        elif r == 2: ret *= -1
        elif r == 3: ret *= 1
        else:        ret *= 0
    return ret


def Tr_Tm_WQ(k, N, m, Q, is_NS):
    assert gcd(Q, N//Q) == 1
    ret = 0
    for N_pm in divisors(N):
        Q_N_pm = gcd(Q, N_pm)
        if is_NS:
            ret += get_alpha(N//N_pm) * s_k_N(k, N_pm, m, Q_N_pm)
        else:
            if not is_squarefree(N // N_pm): continue
            ret += moebius(Q // Q_N_pm) * s_k_N(k, N_pm, m, Q_N_pm) 
    return ret




def test_Tr_Tm_WQ(N_ub, is_NS):
    for N in range(1, N_ub+1):
        print(f'Testing: N={N}')
        for k in range(2, 13, 2):
            print(f'N={N}, k={k}')
            Sk = ModularSymbols(Gamma0(N),k,sign=1).cuspidal_subspace()
            if is_NS: 
                Sk = Sk.new_subspace()
            for Q in divisors(N):
                if gcd(Q, N//Q) != 1: 
                    continue
                for m in range(1, 16):
                    if gcd(m,N) != 1: continue
                    tr = Tr_Tm_WQ(k, N, m, Q, is_NS)
                    WQ = 1/Q^(k//2-1) * Sk.atkin_lehner_operator(Q).matrix()
                    Tm = Sk.hecke_operator(m).matrix()
                    tr_mat = (Tm * WQ).trace() 
                    # print(f'(N={N}, k={k}, Q={Q}) {m}: {tr}, {tr_mat}')
                    assert tr == tr_mat




################################################################################
##### Search sign-pattern weight indexed c2 sequences #############################
################################################################################


def exact_divisors(N):
    return [Q for Q in divisors(N) if gcd(Q, N//Q)==1]
def omega(N):
    return len(factor(N))
def psi(N):
    return product([(p+1)*p^(r-1) for (p,r) in factor(N)])



class signpattern:
    def get_all(N):
        ret = []
        tt = omega(N)
        for mask in range(2^tt):
            sgnpatt_str = ''.join('+-'[(mask>>i)&1] for i in reversed(range(tt)))
            ret.append(signpattern(N, sgnpatt_str))
        return ret
    
    def __init__(self, N, sgnpatt_str):
        sgnpatt_lst = [{'+':1, '-':-1}[c] for c in sgnpatt_str]
        self.d = {p^r: sn for ((p,r),sn) in zip(factor(N), sgnpatt_lst)}
        self.sgnpatt_str = sgnpatt_str
    
    def __repr__(self):
        return self.sgnpatt_str

    def __call__(self, t):
        ret = 1
        for p,r in factor(t):
            ret *= self.d[p^r]
        return ret








    

def Tr_Tm_sgnpatt(k, N, m, sgnpatt, is_NS):
    ret = 0
    for Q in exact_divisors(N):
        ret += sgnpatt(Q) * Tr_Tm_WQ(k, N, m, Q, is_NS=is_NS)
    assert ret % 2^omega(N) == 0
    return ret // 2^omega(N)



# def Tr_Tm_pm_sgnpatt(k, N, m, sgnpatt, is_NS):
#     return m^(-(k-1)/2) * Tr_Tm_sgnpatt(k, N, m, sgnpatt, is_NS)


def c2_Tm_sgnpatt(k, N, m, sgnpatt, is_NS):
    ret = 0
    ret += Tr_Tm_sgnpatt(k, N, m, sgnpatt, is_NS)^2
    for d in divisors(m):
        ret -= d^(k-1) * Tr_Tm_sgnpatt(k, N, (m//d)^2, sgnpatt, is_NS)
    assert ret % 2 == 0
    ret //= 2
    return ret

# def c2_Tm_pm_sgnpatt(k, N, m, sgnpatt, is_NS):
#     ret = 0
#     ret += Tr_Tm_pm_sgnpatt(k, N, m, sgnpatt, is_NS)^2
#     for d in divisors(m):
#         ret -= Tr_Tm_pm_sgnpatt(k, N, (m//d)^2, sgnpatt, is_NS)
#     ret /= 2
#     return ret



def Hecke_alg_mult(A, B, k, N):
    C = {}
    for a in A:
        for b in B:
            assert gcd(a,N) == gcd(b,N) == 1
            for d in divisors(gcd(a,b)):
                n = a*b // d^2
                C[n] = C.get(n, 0) + A[a]*B[b]*d^(k-1)
    return C




def get_Tm_c_coeffs(k, N, m, sgnpatt, is_NS, dim, num_coeffs):
    assert num_coeffs <= dim
    pp = [0]*(num_coeffs+1)
    Tm_r = {1: 1} # Tm^0 = 1 * T_1
    pp[0] = dim   # Tr Tm^0 = dim
    for r in range(1, num_coeffs+1):
        Tm_r = Hecke_alg_mult(Tm_r, {m: 1}, k, N)
        pp[r] = sum(Tm_r[n] * Tr_Tm_sgnpatt(k, N, n, sgnpatt, is_NS) for n in Tm_r)

    cc = [0]*(num_coeffs+1)
    cc[0] = 1
    for r in range(1, num_coeffs+1):
        cc[r] = ZZ(-1/r * sum(cc[r-j]*pp[j] for j in range(1,r+1)))

    if num_coeffs >= 1:
        assert cc[1] == -1 * Tr_Tm_sgnpatt(k, N, m, sgnpatt, is_NS)
    if num_coeffs >= 2:
        assert cc[2] == c2_Tm_sgnpatt(k, N, m, sgnpatt, is_NS)
    return cc



# See Proposition 7.1
def get_k_ub(N, m):
    m_expr = m^5 * sigma(m,0)^3 / sigma(m,1)
    N_expr = psi(N) / (2^omega(N) * N * log(4*N)^2 * sigma(N,0)^4)
    return ceil(((1218 * m_expr) / N_expr) * 24 + 1)

# for N in range(1, 20):
#     print(N, get_k_ub(N,2))
#     print(N, get_k_ub(N,3))
#     print(N, get_k_ub(N,5))
# print('m=2, N=3:', get_k_ub(3,2))







def search_weight_c2_Tm_sgnpatt(m, N, sgnpatt, k_lb, k_ub):
    assert gcd(m,N) == 1
    for k in range(k_lb, k_ub+1, 2):
        if k % 1000 == 2:
            print(f'checking k={k}')
        c2 = c2_Tm_sgnpatt(k, N, m, sgnpatt, is_NS=False)
        if c2 >= 0:
            dim = Tr_Tm_sgnpatt(k, N, 1, sgnpatt, is_NS=False)
            print(f'[dim={dim}] (N,sigma)=({N},{sgnpatt}), k={k}, c2(T{m}) = {c2};{c2/m^(k-1)}')
    print(f'All done with m={m}, (N,sigma)=({N},{sgnpatt})')




def check_conj_general_r(k_ub, N_ub, m_ub, is_NS, max_r=100000):
    COUNTEREXAMPLES = []
    for N in range(1, N_ub+1):
        for sgnpatt in signpattern.get_all(N):
            print(f'CHECKING N={N}, sigma={sgnpatt} ############################')
            DIM = [None] * (k_ub+1)
            for k in range(2, k_ub+1, 2):
                DIM[k] = Tr_Tm_sgnpatt(k, N, 1, sgnpatt, is_NS)
            for m in range(1, m_ub+1):
                if gcd(m, N) != 1:
                    continue
                for k in range(6, k_ub+1, 2):
                    num_coeffs = min(DIM[k], max_r)
                    cc = get_Tm_c_coeffs(k, N, m, sgnpatt, is_NS, DIM[k], num_coeffs)
                    for r in range(1, num_coeffs+1):
                        if cc[r] == 0: 
                            print(f'FOUND: N={N}, {sgnpatt}, m={m}, k={k}, dim={DIM[k]}:  c{r}={cc[r]}')
                            COUNTEREXAMPLES.append((f'c{r}',is_NS,N,sgnpatt,m,k,DIM[k],cc[r]))
    print('Counterexamples:')
    [print(cxmp) for cxmp in COUNTEREXAMPLES]
    print('All done!')



def check_conj_c2(k_ub, N_ub, m_ub, is_NS):
    COUNTEREXAMPLES = []
    for N in range(1, N_ub+1):
        for sgnpatt in signpattern.get_all(N):
            print(f'CHECKING N={N}, sigma={sgnpatt} ############################')
            DIM = [None] * (k_ub+1)
            for k in range(2, k_ub+1, 2):
                DIM[k] = Tr_Tm_sgnpatt(k, N, 1, sgnpatt, is_NS)
            for m in range(1, m_ub+1):
                if gcd(m, N) != 1:
                    continue
                for k in range(8, k_ub+1, 2):
                    if DIM[k] < 2: continue
                    c2 = c2_Tm_sgnpatt(k, N, m, sgnpatt, is_NS)
                    if c2 == 0: 
                        print(f'FOUND: N={N}, {sgnpatt}, m={m}, k={k}, dim={DIM[k]}:  c2={c2}')
                        COUNTEREXAMPLES.append((is_NS,N,sgnpatt,m,k,DIM[k],c2))
    print('Counterexamples:')
    [print(cxmp) for cxmp in COUNTEREXAMPLES]
    print('All done!')



# check_conj_c2(40, 100, 20, True)
# check_conj_c2(100, 1000, 50, False)
# check_conj_general_r(16, 100, 25, True, max_r=3)
# search_weight_c2_Tm_sgnpatt(2, 3, signpattern(3, '-'), 2*10^6, 2*10^6 + 4)
# search_weight_c2_Tm_sgnpatt(2, 15, signpattern(15,'+-'), 5000)
# search_weight_c2_Tm_sgnpatt(5, 6, signpattern(6,'+-'), 5000)












main()
# sage compute_coeffs.sage check_conj_c2 40 300 50 FS
# sage compute_coeffs.sage check_conj_c2 40 300 50 NS
