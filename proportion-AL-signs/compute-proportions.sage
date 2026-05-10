# We implement the trace formula from Skoruppa-Zagier



# This class computes the values of the Lucas Sequence of the first kind.
# This is most efficient when the values for k are accessed sequentially.
class U_seq:
    SEQS = {}
    def seq(t, m):
        if (t, m) not in U_seq.SEQS:
            U_seq.SEQS[ (t, m)  ] = U_seq(t, m)
        return U_seq.SEQS[ (t, m) ]

    def __init__(self, t, m):
        self.t = t
        self.m = m
        self.d = [0,1]
        self.cur = 1

    def get(self,k):
        if k != self.cur: 
            if k < self.cur:
                self.d = [0,1]
                self.cur = 1
            # build up to k
            for i in range(self.cur+1, k+1):
                self.d[i%2] = self.t * self.d[(i-1)%2]  -  self.m * self.d[(i-2)%2]
                self.cur += 1
        # Note, we are asserting the the output is integer (even though t is not)
        return ZZ( self.d[self.cur%2] )
    



def HK_classnum(D):
    return QQ(pari(abs(D)).qfbhclassno())

def H_N(N, Delta):
    a2b = gcd(N, abs(Delta))
    b = squarefree_part(a2b)
    a2 = a2b // b
    if abs(Delta) % (a2*b^2) != 0:
        return 0
    return a2b * kronecker(Delta//(a2*b^2), N//a2b) * HK_classnum(abs(Delta)//(a2*b^2))







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
    ret2 = N // squarefree_part(N)
    return ZZ(sqrt(ret2))
        


def s_k_N(k,N,m,Q):
    term1 = 0
    for Q_pm in divisors(Q):
        for s in get_s(m,Q,Q_pm):
            U_k1 = U_seq.seq(s/sqrt(Q_pm),m).get(k-1) 
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
    assert ret.is_integer()
    return ZZ(ret)



def get_alpha(N):
    ret = 1
    for p,r in factor(N):
        if   r == 1: ret *= -1
        elif r == 2: ret *= -1
        elif r == 3: ret *= 1
        else:        ret *= 0
    return ret

def FS_Tr_Tm_WQ(k, N, m, Q):
    assert gcd(Q, N//Q) == 1
    ret = 0
    for N_pm in divisors(N):
        if is_squarefree(N // N_pm):
            Q_N_pm = gcd(Q, N_pm)
            mob = moebius(Q // Q_N_pm)
            if mob != 0:
               ret += mob * s_k_N(k, N_pm, m, Q_N_pm) 
    return ret


def NS_Tr_Tm_WQ(k, N, m, Q):
    assert gcd(Q, N//Q) == 1
    ret = 0
    for N_pm in divisors(N):
        Q_N_pm = gcd(Q, N_pm)
        ret += get_alpha(N//N_pm) * s_k_N(k, N_pm, m, Q_N_pm)
    return ret




def test_trace(is_NS):
    for N in range(1, 40):
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
                    if is_NS:  tr = NS_Tr_Tm_WQ(k, N, m, Q)
                    else:      tr = FS_Tr_Tm_WQ(k, N, m, Q)
                    WQ = 1/Q^(k//2-1) * Sk.atkin_lehner_operator(Q).matrix()
                    Tm = Sk.hecke_operator(m).matrix()
                    tr_mat = (Tm * WQ).trace() 
                    print(f'(N={N}, k={k}, Q={Q}) {m}: {tr}, {tr_mat}')
                    assert tr == tr_mat



# test_trace(True)



# Now, we compute proportions manually (to compare with predicted proportions)

def print_float_lst(lst):
    str_lst = []
    for v in lst:
        str_lst.append( f'{v.n():.4f}' )
    print( '[' + ', '.join(str_lst) + ']' )

def is_cbfr_square(N):
    return all(r==2 for (p,r) in factor(N))

def prop_formula(N):
    if not is_cbfr_square(N):
        return 1/2

    prod = 1
    for p,r in factor(N):
        prod *= -1 / (p^2 - p - 1)
    return (1 + prod) / 2


def get_one_prop_NS(N, k):
    dim = NS_Tr_Tm_WQ(k, N, 1, 1)
    if dim == 0:
        return -1
    dim_plus = (1/2) * (dim + NS_Tr_Tm_WQ(k, N, 1, N)) 
    return dim_plus / dim



def print_props_NS(N_ub, k_UB):
    for N in range(1, N_ub+1):
        prp = prop_formula(N)
        print(f'N = {N}: prop = {prp.n():.4f}')
        props = [get_one_prop_NS(N,k) for k in range(2, k_UB+1, 2)]
        print_float_lst(props)
        print()


print_props_NS(50, 400)


