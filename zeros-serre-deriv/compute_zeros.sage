import matplotlib.pyplot as plt
import numpy as np
import sys


def main():
    ARG = sys.argv[1]
    if ARG == 'verify_interlacing_Fk_Fsharpk':
        verify_interlacing_Fk_Fsharpk(1000)
    elif ARG == 'verify_loc_Fsharpk_Fhatsharpk':
        verify_loc_Fsharpk_Fhatsharpk(1000)
    elif ARG == 'verify_interlacing_Fk_Fl':
        verify_interlacing_Fk_Fl(1000)
    elif ARG == 'verify_Stieltjes_interlacing_Fk_Fl':
        verify_Stieltjes_interlacing_Fk_Fl(1000)
    elif ARG == 'plot_Fk_minus_Fhatk':
        plot_Fk_minus_Fhatk(10, 20)
    elif ARG == 'plot_B0_B1_B2_B3_U_V':
        plot_B0_B1_B2_B3_U_V()
    elif ARG == 'plot_several_Fk':
        plot_several_Fk(10, 28)
    else:
        assert False, 'incorrect ARG: ' + ARG







pi = RDF.pi()

# Construct Ek, DEk
def qseries_value(F, tau):
    q = CDF(exp(2*pi*I*tau))
    return sum(CDF(F[n]) * q^n for n in range(F.prec()))    
def make_Ek_DEk(k):
    if 2 <= k <= 100:
        PREC = 20 if k == 2 else 100
        Ek_ = -2*k/bernoulli(k) * eisenstein_series_qexp(k,101)
        DEk_ = Ek_.derivative() * Ek_.parent().gen()
        Ek  = lambda z: qseries_value(Ek_, z)
        DEk = lambda z: qseries_value(DEk_, z)
        return Ek, DEk
    else:
        def Ek(z):
            ret = CDF(1)
            for c in range(1,10):
                for d in range(-9, 10):
                    if gcd(c,d) == 1:
                        ret += CDF((c*z+d)^(-k))
            return ret    
        def DEk(z):
            ret = CDF(0)
            for c in range(1,10):
                for d in range(-9, 10):
                    if gcd(c,d) == 1:
                        ret += CDF(c * (c*z+d)^(-k-1))
            return -k / (2*pi*I) * ret
        return Ek, DEk

# Construct B0().
E2, _  = make_Ek_DEk(2)   
def B0(theta):
    z = CDF(exp(I*theta))
    return pi/3 * (z * E2(z)).real() 
B1 = lambda x : sqrt(1 + B0(x)^2)
B2 = lambda x : arctan(B0(x))
B3 = lambda x : -1/2 * (tan(x/2) + B0(x))
dd = lambda x : 2 * cos(x/2)



# Construct Fk and Fhatk
def make_Fk(k, PREC=100):
    Ek, DEk = make_Ek_DEk(k)
    def Fk(theta):
        z = CDF(exp(I*theta))
        ret = DEk(z) - k/12 * E2(z) * Ek(z)
        ret *= z**((k+2)//2)
        return 2*pi/k * ret.real()
    return Fk

def make_Fstark(k):
    def Fstark(x):
        return B1(x) * sin(k/2 * x - B2(x)) 
    return Fstark

def make_Fhatk(k):
    def Fhatk(x):
        return B1(x) * sin(k/2 * x - B2(x)) + B3(x) * dd(x)^(-k)
    return Fhatk

def make_Fsharpk(k, PREC=100):
    Ek2,_ = make_Ek_DEk(k+2)
    def Fsharpk(theta):
        z = CDF(exp(I*theta))
        ret =  Ek2(z)
        ret *= z**((k+2)//2)
        return -1/2 * ret.real()
    return Fsharpk

def make_Fhatsharpk(k):
    def Fhatsharpk(x):
        return -cos((k+2)/2 * x) - 1/2 * dd(x)^(-k-2)
    return Fhatsharpk


def n_k(k):
    n = (k+2) // 12
    if (k+2) % 12 == 2:
        n -= 1
    return n

def delta_k(k):
    if   k % 6 == 0: return 0
    elif k % 6 == 2: return 1
    elif k % 6 == 4: return 2

def presample_zeros(k):
    nk = n_k(k)
    deltak = delta_k(k)
    ret = []
    for i in range(1,nk+1):
        ret.append(2*pi/3 - 2*pi/k * (i - deltak/3))
    return ret

def postsample_zeros(k):
    Fhatk = make_Fhatk(k)
    ret = []
    for t0 in presample_zeros(k):
        try:
            root = find_root(Fhatk, t0 - 1.28 / k, t0 + 0.8 / k)
        except:
            print(f'ERROR: no posts root found; k={k}, t0={t0}, [{t0 - 1.28 / k}, {t0 + 0.8 / k}]')
            continue
        ret.append(root)
    return ret

def actual_zeros(k):
    Fk = make_Fk(k)
    ret = []
    for t0 in presample_zeros(k):
        try:
            root = find_root(Fk, t0 - 1.2 / k, t0 + 1.2 / k)
        except:
            print(f'ERROR: no actual root found; k={k}, t0={t0}, [{t0 - 1.2 / k}, {t0 + 1.2 / k}]')
            continue
        ret.append(root)
    return ret


def sample_zeros(k):
    Fstark = make_Fstark(k)
    ret = []
    for t0 in presample_zeros(k):
        try:
            root = find_root(Fstark, t0 - 1.2 / k, t0 + 1.2 / k)
        except:
            print(f'ERROR: no sample root found; k={k}, t0={t0}, [{t0 - 1.2 / k}, {t0 + 1.2 / k}]')
            continue
        ret.append(root)
    return ret


def actual_zeros_sharp(k):
    Fsharpk = make_Fsharpk(k)
    ret = []
    for t0 in sample_zeros(k):
        try:
            root = find_root(Fsharpk, t0 - 0.25 / k, t0 + 0.25 / k)
        except:
            print(f'ERROR: no actual sharp root found; k={k}, t0={t0}, [{t0 - 0.25 / k}, {t0 + 0.25 / k}]')
            continue
        ret.append(root)
    return ret






###################################################################
#################### Plot functions ###############################
###################################################################

def plot_Fk_Fhatk_presz_postsz(k):
    plt.clf()
    Fk = make_Fk(k)
    Fhatk = make_Fhatk(k)
    for (func, col, nm) in [(Fhatk, 'green', 'Fhat_k'), (Fk, 'purple', 'F_k')]:
        xx = np.linspace(np.pi/2, 2*np.pi/3, 10*k) 
        yy = [func(xx_) for xx_ in xx]
        plt.plot(xx, yy, color=col, label=nm)
    for (xx, col, nm) in [(presample_zeros(k), 'blue', 'presz'), (postsample_zeros(k), 'red', 'postsz')]:
        yy = [0] * len(xx)
        plt.scatter(xx, yy, color=col, marker='x', s=20)
    plt.legend(loc='lower center')
    plt.grid(True)
    plt.show()


def plot_Fhatsharpk_Fhatk_presz_postsz(k):
    plt.clf()
    Fhatk = make_Fhatk(k)
    Fhatsharpk2 = make_Fhatsharpk(k)
    for (func, col, nm) in [(Fhatsharpk2, 'orange', 'Fhatsharp_k'), (Fhatk, 'green', 'Fhat_k')]:
        xx = np.linspace(np.pi/2, 2*np.pi/3, 10*k) 
        yy = [func(xx_) for xx_ in xx]
        plt.plot(xx, yy, color=col, label=nm)
    for (xx, col, nm) in [(presample_zeros(k), 'blue', 'presz'), (postsample_zeros(k), 'red', 'postsz')]:
        yy = [0] * len(xx)
        plt.scatter(xx, yy, color=col, marker='x', s=20)
    plt.legend(loc='lower center')
    plt.grid(True)
    plt.show()


def plot_Fhatsharpk_minus_Fhatk(k):
    plt.clf()
    Fhatk = make_Fhatk(k)
    Fhatsharpk2 = make_Fhatsharpk(k)
    diff = lambda x: Fhatsharpk2(x) - Fhatk(x) 
    for (func, col, nm) in [(diff, 'orange', 'Fhatsharp_k+2 - Fhat_k')]:
        xx = np.linspace(np.pi/2, 2*np.pi/3, 10*k) 
        yy = [func(xx_) for xx_ in xx]
        plt.plot(xx, yy, color=col, label=nm)
    plt.legend(loc='lower center')
    plt.grid(True)
    plt.show()


def plot_Fk_minus_Fhatk(k_lb, k_ub):
    plt.clf()
    for k in range(k_lb, k_ub+1, 2):
        Fk = make_Fk(k)
        Fhatk = make_Fhatk(k)
        xx = np.linspace(np.pi/2, 2*np.pi/3, 100) 
        yy = [Fk(xx_)-Fhatk(xx_) for xx_ in xx]
        plt.plot(xx, yy, label=f'F{k} - Fhat{k}')
    plt.legend(loc='lower center')
    plt.grid(True)
    plt.show()




def plot_several_Fk(k_lb, k_ub):
    plt.clf()
    for k in range(k_lb, k_ub+1, 2):
        Fk = make_Fk(k)
        for (func, nm) in [(Fk, f'F{k}')]:
            xx = np.linspace(np.pi/2, 2*np.pi/3, 20*k) 
            yy = [func(xx_) for xx_ in xx]
            plt.plot(xx, yy, label=nm)
    plt.legend(loc='lower center')
    plt.grid(True)
    plt.show()



def plot_B0_B1_B2_B3_U_V():
    def deriv(func, x, dx=1e-5):
        return (func(x+dx) - func(x-dx)) / (2*dx)
    plt.clf()

    # 0th deriv
    # B0, B1, B2, B3 already defined
    UU = lambda x : 2 * B2(x)
    VV = lambda x : UU(x) * deriv(UU, x)
    for (func, nm) in [(B0,'B0'), (B1,'B1'), (B2,'B2'), (B3,'B3'), (UU,'U'), (VV,'V')]:
        xx = np.linspace(np.pi/2, 2*np.pi/3, 1000) 
        yy = [func(xx_) for xx_ in xx]
        plt.plot(xx, yy, label=nm)
    plt.legend(loc='lower center')
    plt.grid(True)
    plt.show()

    # 1st deriv
    B0pm = lambda x : deriv(B0, x)
    B1pm = lambda x : deriv(B1, x)
    B2pm = lambda x : deriv(B2, x)
    B3pm = lambda x : deriv(B3, x)
    UUpm = lambda x : deriv(UU, x)
    VVpm = lambda x : deriv(VV, x)
    for (func, nm) in [(B0pm,"B0'"), (B1pm,"B1'"), (B2pm,"B2'"), (B3pm,"B3'"), (UUpm,"U'"), (VVpm,"V'")]:
        xx = np.linspace(np.pi/2, 2*np.pi/3, 1000) 
        yy = [func(xx_) for xx_ in xx]
        plt.plot(xx, yy, label=nm)
    plt.legend(loc='lower center')
    plt.grid(True)
    plt.show()

    # 2nd deriv
    B0pmpm = lambda x : deriv(B0pm, x)
    B2pmpm = lambda x : deriv(B2pm, x)
    UUpmpm = lambda x : deriv(UUpm, x)
    for (func, nm) in [(B0pmpm,"B0''"), (B2pmpm,"B2''"), (UUpmpm,"U''")]:
        xx = np.linspace(np.pi/2, 2*np.pi/3, 400) 
        yy = [func(xx_) for xx_ in xx]
        plt.plot(xx, yy, label=nm)
    plt.legend(loc='lower center')
    plt.grid(True)
    plt.show()
    




############################################################################
####################### verify Interlacing claims    #######################
############################################################################




def check_Stieltjes_interlacing(pts1, pts2, EPS):
    if len(pts2) <= 1:
        return True
    if len(pts1) == 0 or pts1[-1] >= pts2[-2]-EPS:
        return False
    
    idx = 0
    for i in range(len(pts2)-1):
        ub = pts2[i]
        lb = pts2[i+1]
        # increment idx until pts1[idx] < ub-EPS
        while pts1[idx] >= ub-EPS: idx += 1
        if not (lb+EPS < pts1[idx] < ub-EPS):
            return False
    return True


def check_interlacing(pts1, pts2, EPS):
    return check_Stieltjes_interlacing(pts1, pts2, EPS) and check_Stieltjes_interlacing(pts2, pts1, EPS)






EXCEPTIONS = [
            (4, 14), (4, 18), (4, 20), (4, 24), (6, 16), (6, 20), (6, 24), (8, 18), (8, 24), 
            (10, 20), (10, 24), (10, 26), (10, 28), (10, 30), (10, 32), (10, 36), (14, 24), (14, 30), 
            (14, 32), (14, 36), (16, 32), (20, 30), (20, 36), (22, 38), (26, 42), (32, 48), (38, 54)
]
def verify_interlacing_Fk_Fl(k_ub):
    ZEROS = [[] for i in range(k_ub+27)]
    for k in range(4, k_ub+27, 2):
        if k % 100 == 4:
            print(f'Computing zeros for k={k}...')
        ZEROS[k] = actual_zeros(k)
    for k in range(4, k_ub+1, 2):
        for ell in range(k+2, k+27, 2):
            cond1 = bool(ell-k in [2,4,6,8,12]  or  (k,ell) in EXCEPTIONS)
            cond2 = check_interlacing(ZEROS[ell], ZEROS[k], 0)
            assert cond1 == cond2
    print(f'Done with verifying interlacing up to {k_ub}.')




def verify_Stieltjes_interlacing_Fk_Fl(k_ub):
    ell_ub = lambda k: round(2.3*k)+2
    ZEROS = [[] for i in range(ell_ub(k_ub)+1)]
    for k in range(4, ell_ub(k_ub)+1, 2):
        if k % 100 == 4:
            print(f'Computing zeros for k={k}...')
        ZEROS[k] = actual_zeros(k)
    for k in range(4, k_ub+1, 2):
        for ell in range(k+2, ell_ub(k)+1, 2):
            assert check_Stieltjes_interlacing(ZEROS[ell], ZEROS[k], 10^(-11))
    print(f'Done with verifying Stieltjes interlacing up to {k_ub}.')




def verify_interlacing_Fk_Fsharpk(k_ub):
    for k in range(4, k_ub+1, 2):
        zeros1 = actual_zeros(k)
        zeros2 = actual_zeros_sharp(k)
        if k % 100 == 4:
            print(f'verifying Fk [{len(zeros1)}], Fsharpk [{len(zeros2)}] zeros for k={k}...')
        assert check_interlacing(zeros1, zeros2, 10^(-8))
    print(f'Done with verifying interlacing Fk Fsharpk up to {k_ub}.')





def verify_loc_Fsharpk_Fhatsharpk(k_ub):
    FAILURE = False
    EPS = 10^(-11)
    for k in range(4, k_ub+1, 2):
        zeros1 = sample_zeros(k)
        if k % 100 == 4:
            print(f'verifying Fsharpk,Fhatsharpk [{len(zeros1)}] zeros for k={k}...')
        for t0 in zeros1:
            for func in [make_Fsharpk(k), make_Fhatsharpk(k)]:
                # Existence -------------
                try:
                    root = find_root(func, t0 - 0.25 / k, t0 + 0.25 / k)
                except:
                    print(f'ERROR: no sharp root found; k={k}, t0={t0}, [{t0 - 0.25 / k}, {t0 + 0.25 / k}]')
                    FAILURE = True
                # Uniqueness ---------------
                try:
                    root_left = find_root(func, t0 - 0.25 / k, root-EPS)
                    print(f'ERROR: sharp root_left found; k={k}, t0={t0}, [{t0 - 0.25 / k}, {root-EPS}]: {root_left}')
                    FAILURE = True
                except:
                    pass
                try:
                    root_right = find_root(func, root+EPS, t0 + 0.25 / k)
                    print(f'ERROR: sharp root_right found; k={k}, t0={t0}, [{root+EPS}, {t0+0.25/k}]: {root_right}')
                    FAILURE = True
                except:
                    pass
                assert not FAILURE
    print(f'Done with verifying location up to {k_ub}.')

                



main()



