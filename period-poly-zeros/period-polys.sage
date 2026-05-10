from sage.modular.dims import dimension_new_cusp_forms
import sys



X = PolynomialRing(CC, 'X').gen()

def load_L_values(filename):
    ret = {}
    with open(filename, 'r') as f:
        lns = f.readlines()
    idx = 0
    while idx < len(lns):
        Nkdim = eval(lns[idx]); idx += 1
        N, k, dim = [ZZ(int(Nkdim[var_str])) for var_str in ('N','k','dim')]
        assert dim == dimension_new_cusp_forms(N,k)
        data = []
        for i in range(dim):
            arr_str = lns[idx]; idx += 1
            L_vals = eval(arr_str.replace(' E','E'))
            assert type(L_vals) == list
            assert len(L_vals) == k-1
            data.append([None] + L_vals)
        ret[(N,k)] = data
    return ret
    
L_VALUE_DATA = load_L_values(sys.argv[1])


def per_poly(L_vals, k, plus):
    if plus: n_inds = range(0, (k-2)+1, 2)
    else:    n_inds = range(1, (k-3)+1, 2)

    ret = 0
    for n in n_inds:
        ret += (2*pi.n()*i * X)^n / factorial(n) * L_vals[k-n-1]
    return ret




def get_nonzero_roots(poly, EPS):
    deg = poly.degree()
    zero_cnt = 0
    ret = []
    for rt,mult in poly.roots(CC):
        if abs(rt) < EPS:
            zero_cnt += mult
        else:
            for i in range(mult):
                ret.append(rt)

    if deg % 2 == 1:
        assert zero_cnt == 1
        assert len(ret) == deg-1
    else:
        assert zero_cnt == 0
        assert len(ret) == deg
    return ret


def list_diff(lst1, lst2): 
    if len(lst1) != len(lst2):
        return float('inf')
    lst1.sort()
    lst2.sort()
    diff = [abs(lst1[i]-lst2[i]) for i in range(len(lst1))]
    return max(diff)




def beta(N, plus):
    if plus:
        if   N == 1:       return 2
        elif 2 <= N <= 15: return 1
        else:              return 0
    else:
        if   1 <= N <= 3:  return 1
        else:              return 0



def check_zero_locations(N, k, plus, condition, loc_vals=None):
    EPS = 10^(-7)
    for j,L_vals in enumerate(L_VALUE_DATA[(N,k)]):
        ply = per_poly(L_vals, k, plus)
        roots = get_nonzero_roots(ply, EPS)
        # print(roots)
        roots_off_circle = []
        for rt in roots:
            dst = abs(rt)
            if abs(dst - 1/sqrt(N)) >= EPS:
                roots_off_circle.append(rt)
        
        if condition == 'print diff loc_vals':
            diff = list_diff(roots_off_circle, loc_vals)
            print(diff)
            if diff == float('inf'):
                print(roots_off_circle)
            continue
        
        if condition == '>4beta':
            cond_sat = ( len(roots_off_circle) > 4*beta(N, plus) )
        elif condition == '==4beta':
            cond_sat = ( len(roots_off_circle) == 4*beta(N, plus) )
        elif condition == '1/4 not a root':
            cond_sat = min(abs(rt-1/4) for rt in roots) > EPS
        else: 
            assert False        
        
        if cond_sat:
            print(f'FOUND EXAMPLE: ({condition})   N={N}, k={k}, plus={plus}, j={j}')
            print(roots_off_circle)
            print(f'1/sqrt(N) = {1/sqrt(N).n()}')
            print(ply)
            print()







########################################################
########################################################



def k_N(N, plus):
    if plus:
        if   N == 1:         return 76
        elif N == 2:         return 60
        elif N == 3:         return 52
        elif  4 <= N <= 5:   return 44
        elif  6 <= N <= 18:  return 36
        elif 19 <= N <= 29:  return 28
        elif 30 <= N <= 103: return 20
        else: assert False
    else:
        if    2 <= N <= 4:    return 56
        elif  5 <= N <= 7:    return 40
        elif  8 <= N <= 13:   return 32
        elif 14 <= N <= 39:   return 24
        elif 40 <= N <= 2499: return 16
        else: assert False



def verify_theorem_plus():
    for N in range(1,104):
        for k in range(12, k_N(N,True), 2):
            if N == 5: continue
            print(f'thm [plus]  N={N} \t  k={k}:')
            check_zero_locations(N, k, True, '>4beta')
    
def verify_theorem_minus():
    for N in range(2,2500):
        for k in range(2*ceil(log(N,2)), k_N(N,False), 2):
            print(f'thm [minus]  N={N} \t  k={k}:')
            check_zero_locations(N, k, False, '>4beta')

def verify_optimality_kN_plus():
    for N in [17]:
        for k in range(6, 12, 2):
            print(f'optk [plus]  N={N} \t  k={k}:')
            check_zero_locations(N, k, True, '>4beta')

def verify_optimality_kN_minus():
    for k in range(6, 50, 2):
        for N in range(1, 130):
            if (N,k) not in L_VALUE_DATA: continue
            print(f'optk [minus]  N={N} \t  k={k}:')
            check_zero_locations(N, k, True, '>4beta')

def verify_optimality_betaN_plus():
    for k in [24]:
        for N in range(1,16):
            print(f'optb [plus]  N={N} \t  k={k}:')
            check_zero_locations(N, k, True, '==4beta')

def verify_optimality_betaN_minus():
    for k in [24]:
        for N in range(1,4):
            print(f'optb [minus]  N={N} \t  k={k}:')
            check_zero_locations(N, k, False, '==4beta')

def verify_location_values_plus():
    for N in range(1, 16):
        loc_vals = [1/4, -1/4, 4/N, -4/N]
        if N == 1: 
            loc_vals.extend([3/4,-3/4,4/3,-4/3])
        for k in range(6, 102, 2):
            if (N,k) not in L_VALUE_DATA: continue
            print(f'locv [plus]  N={N} \t  k={k}:')
            check_zero_locations(N, k, True, 'print diff loc_vals', loc_vals=loc_vals)

def verify_location_values_minus():
    for N in range(1, 4):
        loc_vals = [1/2, -1/2, 2/N, -2/N]
        for k in range(6, 102, 2):
            if (N,k) not in L_VALUE_DATA: continue
            print(f'locv [minus]  N={N} \t  k={k}:')
            check_zero_locations(N, k, False, 'print diff loc_vals', loc_vals=loc_vals)

def verify_conjecture():
    for N in [16]:
        for k in range(4, 102, 2):
            print(f'conj      N={N} \t  k={k}:')
            check_zero_locations(N, k, True, '1/4 not a root')


ARG = sys.argv[2]
if   ARG == 'theorem_plus':
    verify_theorem_plus()
elif ARG == 'theorem_minus':
    verify_theorem_minus()
elif ARG == 'optimality_kN_plus':
    verify_optimality_kN_plus()
elif ARG == 'optimality_kN_minus':
    verify_optimality_kN_minus()
elif ARG == 'optimality_betaN_plus':
    verify_optimality_betaN_plus()
elif ARG == 'optimality_betaN_minus':
    verify_optimality_betaN_minus()
elif ARG == 'location_values_plus':
    verify_location_values_plus()
elif ARG == 'location_values_minus':
    verify_location_values_minus()
elif ARG == 'conjecture':
    verify_conjecture()
else:
    assert False


