import sys
from math import sqrt

K = 584000001


# sieve[N] == the largest prime divisor of N
sieve = [0]*K
for p in range(2,K):
    if sieve[p] == 0: # prime 
        for i in range(p, K, p):
            sieve[i] = p


print('done with sieve')
sys.stdout.flush()

def factor(N):
    ret = []
    while N != 1:
        pm = sieve[N]
        exp = 0
        while sieve[N] == pm:
            N //= pm
            exp += 1
        ret.append((pm, exp))
    ret.reverse()
    return ret

def _2_omega(N):
    return 2 ** len(factor(N))

def psi(N):
    ret = 1
    for pm, exp in factor(N):
        ret *= (pm+1) * pm**(exp-1)
    return ret

def theta1(N):
    return sqrt(N) * _2_omega(N) / psi(N)

def theta2(N):
    return _2_omega(N)**2 / psi(N)

def theta3(N):
    return _2_omega(N) / psi(N)


theta1_bds = [100.0, 1.0, 0.465, 0.257, 0.133, 0.0607, 0.0265, 0.0106]
theta2_bds = [100.0, 1.34, 0.445, 0.149, 0.0424, 0.00941, 0.00189, 0.000314]
theta3_bds = [100.0, 1.00, 0.0556, 0.00926, 0.00133, 0.000147, 0.000015, 0.000015]



idx1 = len(theta1_bds) - 1
idx2 = idx1
idx3 = idx1


for N in range(K-1, 0, -1):
    while theta1(N) > theta1_bds[idx1]:
        print(f'N min for: ({idx1}) theta1 <= {theta1_bds[idx1]:.7f}:     {N+1}')
        idx1 -= 1
    while theta2(N) > theta2_bds[idx2]:
        print(f'N min for: ({idx2}) theta2 <= {theta2_bds[idx2]:.7f}:     {N+1}')
        idx2 -= 1
    while theta3(N) > theta3_bds[idx3]:
        print(f'N min for: ({idx3}) theta3 <= {theta3_bds[idx3]:.7f}:     {N+1}')
        idx3 -= 1
    



    


