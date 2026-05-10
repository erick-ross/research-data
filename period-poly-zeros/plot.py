import sys
import matplotlib.pyplot as plt
from math import pi, exp, log, sqrt, factorial
from math import cos, sin, tan, tanh, cosh, sinh
from math import floor, ceil
from math import atan as arctan, acos as arccos
from cmath import exp as cexp, cos as ccos, sin as csin, phase as carg


def cot(x):
    return tan(pi/2-x)



def get_r(theta, N, plus):
    phi = (2*pi) / sqrt(N)
    C0 = phi * cos(theta)
    S0 = phi * sin(theta)
    if plus:
        return sqrt( (1/2)*cosh(2*S0) + (1/2)*cos(2*C0) )
    else:
        return sqrt( (1/2)*cosh(2*S0) - (1/2)*cos(2*C0) )



def get_a_hat(theta, N, plus):
    phi = (2*pi) / sqrt(N)
    C0 = phi * cos(theta)
    S0 = phi * sin(theta)
    if plus:
        return arctan( -tan(C0) * tanh(S0) )
    else:
        return arctan(  cot(C0) * tanh(S0) )




# 0.0 | 3.2 | 8.8 | 10.0
#     0     1     2
def range_idx(x, endpts):
    for i in range(len(endpts)-1):
        if endpts[i] <= x < endpts[i+1]:
            return i
    if x == endpts[-1]:
        return len(endpts)-2
    assert False, 'out of range'


def get_a(theta, N, plus):
    if plus:
        if N >= 17 or N == 16:
            endpts = [0, pi]
            pi_mult = [0]
        elif 2 <= N <= 15:
            endpts = [0, arccos(sqrt(N)/4), arccos(-sqrt(N)/4), pi]
            pi_mult = [1, 2, 3]
        elif N == 1:
            endpts = [0, arccos(3/4), arccos(1/4), arccos(-1/4), arccos(-3/4), pi]
            pi_mult = [0, 1, 2, 3, 4]
    else:
        if N >= 5 or N == 4:
            endpts = [0, pi/2, pi]
            pi_mult = [0, 1]
        elif 2 <= N <= 3:
            endpts = [0, arccos(sqrt(N)/2), pi/2, arccos(-sqrt(N)/2), pi]
            pi_mult = [1, 2, 3, 4]
    a_hat = get_a_hat(theta, N, plus)
    return a_hat + pi_mult[range_idx(theta, endpts)] * pi


def get_a_endpt_vals(N, plus):
    if plus:
        if N >= 17:
            return [0, 0]
        elif  N == 16:
            return [-pi/2, pi/2]
        elif 2 <= N <= 15:
            return [pi, 3*pi]
        elif N == 1:
            return [0, 4*pi]
    else:
        if N >= 5:
            return [0, pi]
        elif N == 4:
            return [-pi/2, 3*pi/2]
        elif 2 <= N <= 3:
            return [pi, 4*pi]

        





def get_S_NMk(N, M, k, plus, verbose=False):
    assert k % 2 == 0
    m = (k-2) // 2
    phi = 2*pi/sqrt(N)

    sum1, sum2, sum3 = 0.0, 0.0, 0.0
    if plus:
        for j in range(0, floor((m-1)/4)+1):
            sum1 += phi**(2*j) / factorial(2*j)
        for j in range(ceil(m/4), floor(m/2)+1):
            weight = 1 if 2*j < m else 1/2 
            sum2 += weight * phi**(2*j) / factorial(2*j)
        sum3 = cosh(phi) - sum1
    else:
        for j in range(0, floor((m-3)/4)+1):
            sum1 += phi**(2*j+1) / factorial(2*j+1)
        for j in range(ceil((m-2)/4), floor((m-1)/2)+1):
            weight = 1 if 2*j+1 < m else 1/2 
            sum2 += weight * phi**(2*j+1) / factorial(2*j+1)
        sum3 = sinh(phi) - sum1
    
    # S1 = zeta2_minus_1((k+2)/4) * sum1
    S1 = 2.71 * 2**(-k/4) * sum1
    S2 = 4 / (exp(pi/sqrt(N)) - 1) * 2**(k/2)
    S2 *= exp(-pi*k/sqrt(M)) # the only place where M appears
    S2 += 4*sqrt(k)*(1 + log(k))
    S2 *= sum2
    S3 = sum3

    if verbose:
        print(f'[N={N},M={M}]     k={k}:')
        print(f'sum123:    {sum1}\t\t{sum2}\t\t{sum3}')
        print(f'S123:      {S1  }\t\t{S2  }\t\t{S3  }')
        print()
    return S1 + S2 + S3



##########################################################################################
##########################################################################################

def arg_correct(val, theta, N, plus):
    if plus:
        true_val = carg(ccos(2*pi*cexp(1j*theta)/sqrt(N)))
    else:
        true_val = carg(csin(2*pi*cexp(1j*theta)/sqrt(N)))
    diff = val - true_val
    while diff < -pi: diff += 2*pi
    while diff >  pi: diff -= 2*pi
    return abs(diff) < 10**(-8)


def rad_correct(val, theta, N, plus):
    if plus:
        true_val = abs(ccos(2*pi*cexp(1j*theta)/sqrt(N)))
    else:
        true_val = abs(csin(2*pi*cexp(1j*theta)/sqrt(N)))
    diff = val - true_val
    return abs(diff) < 10**(-8)


##########################################################################################
##########################################################################################


def plot(title, x_vals, y_vals, y_vals2=None, endpts=None):
    plt.clf()
    if endpts != None:
        plt.plot(endpts[0], endpts[1], 'x', color='green', markersize=10)
    if y_vals2 != None:
        plt.plot(x_vals, y_vals2, '.', color='red', markersize=3)
    plt.plot(x_vals, y_vals)
    plt.title(title)
    plt.grid(True) 
    plt.show()



def interval_mesh(aa, bb, nn):
    delta = (bb-aa)/nn
    return [(aa + i*delta) for i in range(nn+1)]
MESH_0_PI = interval_mesh(0 + 10**(-6), pi - 10**(-6), 1000) 



def plot_arg(plus):
    x_vals = MESH_0_PI
    N_vals = list(range(1, 21)) + [25, 40, 100, 1000, 100000]
    if not plus: N_vals.remove(1)
    
    for N in N_vals:
        y_vals = [get_a(theta, N, plus) for theta in x_vals]
        assert all(arg_correct(val, theta, N, plus) for (theta,val) in zip(x_vals, y_vals))
        y_vals2 = [get_a_hat(theta, N, plus) for theta in x_vals]
        endpts = ([0,pi], get_a_endpt_vals(N,plus))
        title = f'a^{"-+"[plus]}(theta) for N={N}'
        plot(title, x_vals, y_vals, y_vals2=y_vals2, endpts=endpts)


def plot_rad(plus):
    x_vals = MESH_0_PI
    N_vals = list(range(1, 21)) + [25, 40, 100, 1000, 100000]
    if not plus: N_vals.remove(1)

    for N in N_vals:
        y_vals = [get_r(theta, N, plus) for theta in x_vals]
        assert all(rad_correct(val, theta, N, plus) for (theta,val) in zip(x_vals, y_vals))
        title = f'r^{"-+"[plus]}(theta) for N={N}'
        plot(title, x_vals, y_vals)



def plot_SNMk(plus):
    if plus:
        N_M_klb = [
            (1,1,76), (2,2,60), (3,3,52), (4,4,44), (5,5,44), (6,6,36), (7,7,36), (8,8,36), 
            (9,9,36), (10,10,36), (11,11,36), (12,12,36), (13,13,36), (14,14,36), (15,15,36), 
            (17,18,36), (19,29,28), (30,39,20),(40,103,20),
            (104,104,12), (105,109,12),
            (110,119,12),(120,169,12),(170,249,12),
            (250,399,12), (400,699,12), (700,1599,12), (1600,float('inf'),12)
        ]
    else:
        N_M_klb = [
            (2,2,56), (3,3,56), (5,5,40), (6,6,40), (7,7,40), (8,8,32), (9,9,32), (10,10,32), 
            (11,11,32), (12,12,32), (13,13,32), (14,14,24), (15,15,24), 
            (16, 39, 24), 
            (40,49,16),  (50,109,16), (110,279,16), (280,699,16), (700, 2499, 16)
        ]
        

    for N,M,klb in N_M_klb:
        x_vals = list(range(klb, klb+100, 2))
        y_vals = [get_S_NMk(N, M, k, plus) for k in x_vals]
        title = f'S^{"-+"[plus]}_NMk for N={N},M={M},  klb={klb}\n'
        if plus: 
            y_vals2 = [abs(cos(2*pi/sqrt(N))) for k in x_vals]
            title += f'and |cos(2pi/sqrt({N}))| = {y_vals2[0]:.5f}'
        else:    
            y_vals2 = [abs(sin(2*pi/sqrt(M))) for k in x_vals]
            title += f'and |sin(2pi/sqrt({M}))| = {y_vals2[0]:.5f}'
        plot(title, x_vals, y_vals, y_vals2=y_vals2)


def plot_mSNNk(plus):
    if plus:
        klb = 36
        N0 = 16
    else:
        klb = 56
        N0 = 4

    x_vals = list(range(klb, klb+100, 2))
    y_vals = [((k-2)//2) * get_S_NMk(N0, N0, k, plus) for k in x_vals]
    y_vals2 = [1/2 for k in x_vals]
    title = f'm*S^{"-+"[plus]}_NNk for N={N0},  klb={klb}\n'
    plot(title, x_vals, y_vals, y_vals2=y_vals2)
    


##########################################################################################
##########################################################################################

PARAM = sys.argv[1]
PARAM2 = sys.argv[2]
is_plus = {'plus': True, 'minus': False}[PARAM2]
if   PARAM == 'arg':
    plot_arg(is_plus)
elif PARAM == 'rad':
    plot_rad(is_plus)
elif PARAM == 'SNMk':
    plot_SNMk(is_plus)
elif PARAM == 'mSNNk':
    plot_mSNNk(is_plus)
else:
    assert False, 'invalid parameter: ' + PARAM

