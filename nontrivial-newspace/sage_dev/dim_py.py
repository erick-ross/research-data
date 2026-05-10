from sage.all_cmdline import *

from sage.arith.misc import factor
from sage.modular.dirichlet import DirichletGroup
from sage.arith.misc import divisors
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.arith.misc import kronecker
from sage.arith.misc import CRT
from sage.arith.misc import moebius
from sage.rings.integer_ring import ZZ
from sage.modular.dims import dimension_new_cusp_forms


class MyClass:

    # An internal function to implement the explicit dimension formula
    # given in Theorem 1.4 of https://arxiv.org/abs/2407.08881
    def _dimension_new_cusp_forms_Ross(N,k,chi):
        # First term of explicit dimension formula
        def psi_local(p,r):
            assert r >= 1
            return (p+1)*p**(r-1)
        
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
                    return p**2-p-1
                else:
                    return (p**3-p**2-p+1) * p**(r-3)
            else:
                if r==1:
                    return p-2
                else:
                    return (p**2-2*p+1) * p**(r-2)
                
        def beta_psi_f(n, f):
            ret = 1
            for p,r in factor(n):
                ret *= beta_psi_f_local(p,r,f.valuation(p))
            return ret


        # Second term of explicit dimension formula
        def beta_sigma_f_local(p,r,alpha):
            assert r>=1
            if alpha == 0:
                if r % 2 == 1:
                    return 0
                elif r == 2:
                    return p-2
                else:
                    return (p**2-2*p+1) * p**(r//2-2)
            elif r == 1:
                if alpha == 1:
                    return (p-3)/2
                else:
                    return p-2
            else:
                if r >= alpha+1 and (r+alpha)%2==1:
                    return 0
                elif r >= alpha+2 and (r+alpha)%2==0:
                    return ZZ(1)/2 * (p**2-2*p+1) * p**((r+alpha)/2-2) 
                elif r == alpha:
                    return ZZ(1)/2 * (p**2-3*p+2) * p**(r-2)
                else:
                    return (p**2-2*p+1) * p**(r-2) 

        def beta_sigma_f(n,f):
            ret = 1
            for p,r in factor(n):
                ret *= beta_sigma_f_local(p,r,f.valuation(p))
            return ret

        
        # Third term of explicit dimension formula
        def get_chi_p_alpha(f_fact, chi, p, x):
            rems = [(int(x) if q == p else 1) for (q,alpha) in f_fact]
            mods = [p**alpha for (p,alpha) in f_fact]
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
                u = Zmod(p**r)(-3).sqrt(extend=False)
                chi_x = get_chi_p_alpha(f_fact, chi, p, (-1+u)/2)
                assert chi_x**3 == 1
                if chi_x == 1:
                    return 2
                else:
                    return -1
        
        def rho(n,chi,f):
            f_fact = factor(f)
            ret = 1
            for p,r in factor(n):
                ret *= rho_local(p,r,f.valuation(p),chi,f_fact)
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
                ret *= beta_rho_f_local(p,r,f.valuation(p))
            return ret



        # Fourth term ##################################
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
                upm = Zmod(p**r)(-1).sqrt(extend=False)
                chi_x = get_chi_p_alpha(f_fact, chi, p, upm)
                assert chi_x**4 == 1
                if chi_x == 1:
                    return 2
                elif chi_x == -1:
                    return -2
                else:
                    return 0
                
        def rhopm(n,chi,f):
            f_fact = factor(f)
            ret = 1
            for p,r in factor(n):
                ret *= rhopm_local(p,r,f.valuation(p),chi,f_fact)
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
                ret *= beta_rhopm_f_local(p,r,f.valuation(p))
            return ret



        # Constant factors for each term of the explicit dimension formula
        def c3(k):
            return (k-1)/3 - k//3

        def c4(k):
            return (k-1)/4 - k//4

        def c0(k,f):
            if k == 2 and f == 1:
                return 1
            else:
                return 0

        def omega(n):
            return len(factor(n))

            

        # Finally, compute the actual dimension formula
        assert k >= 2
        if chi(-1) != (-1)**k:
            return 0
        f = chi.conductor()
        ret = 0
        ret += (k-1)/12 * psi(f) * beta_psi_f(N//f, f)
        ret += -1 * c3(k) * rho(f, chi, f) * beta_rho_f(N//f, f)
        ret += -1 * c4(k) * rhopm(f, chi, f) * beta_rhopm_f(N//f, f)
        ret += ZZ(-1)/2 * 2**omega(f) * beta_sigma_f(N//f, f) 
        ret += c0(k,f) * moebius(N//f)
        assert ret.is_integral()
        return ZZ(ret)


################################################################


def test_dimension_new_cusp_forms_Ross(N_ub=50, k_ub=10):
    for N in range(1,N_ub+1):
        print('testing', N)
        for eps_idx, eps in enumerate(DirichletGroup(N)):
            for k in range(2,k_ub+1):
                Ross_dim = MyClass._dimension_new_cusp_forms_Ross(ZZ(N),ZZ(k),eps)
                CO_dim = dimension_new_cusp_forms(eps, k)
                assert Ross_dim == CO_dim





import sys
if sys.argv[1] == 'test':
    test_dimension_new_cusp_forms_Ross(int(sys.argv[2]))

