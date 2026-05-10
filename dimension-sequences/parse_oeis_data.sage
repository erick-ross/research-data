import sys
DATA = []
for ln in sys.stdin.readlines():
    DATA.append(int(ln.split()[1]))



def delta_M(M):
    if (-M) % 4 == 1: return -M
    else:             return -4*M



def get_h_pm(delta):
    assert delta < 0 and delta % 4 in [0,1]
    if delta == -3: return 1/3
    if delta == -4: return 1/2

    discr_idx = abs(delta) // 2 - 1
    return DATA[discr_idx]


N_vals = [1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 21, 22, 26, 30, 33, 
            34, 35, 38, 39, 42, 46, 66, 70, 78, 102, 105, 110, 114, 
            130, 138, 210, 330, 390] + [29, 31, 37]


h_pm_dic = {}
for N in N_vals:
    for M in divisors(N):
        delta = delta_M(M)
        h_pm_dic[delta] = get_h_pm(delta)
print(h_pm_dic)


