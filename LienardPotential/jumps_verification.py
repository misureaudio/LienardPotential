import numpy as np
from scipy.integrate import solve_ivp
for MU in [50.0, 100.0]:
    def vdp(t, z):
        x, y = z
        return [y, MU * (1 - x*x) * y - x]
    T_END = 2*162.0
    sol = solve_ivp(vdp, [0.0, T_END], [1.0, 0.0], method='Radau', rtol=1e-10, atol=1e-10, dense_output=True, max_step=0.05)
    t = np.linspace(0.0, T_END, 1200000)
    x = sol.sol(t)[0]; y = sol.sol(t)[1]
    dy = MU * (1 - x*x) * y - x
    jump = np.abs(dy) > 1e3
    idx = np.flatnonzero(jump)
    segs = []
    start = idx[0]
    for a, b in zip(idx[:-1], idx[1:]):
        if b - a > 1:
            segs.append((start, a)); start = b
    segs.append((start, idx[-1]))
    Ce = 2*MU/3
    print(f'mu={MU}: expected |C| = 2mu/3 = {Ce:.3f};  max|x| = {x.max():.4f}, min|x| = {x.min():.4f}')
    for (a, b) in segs:
        if b - a < 5: continue
        s = slice(a+2, b-1)
        C = y[s] - MU*(x[s] - x[s]**3/3)
        # which fold does this jump start from? sign of x at the start
        print(f'  jump x: {x[a]:7.3f} -> {x[b]:7.3f}   C_mean = {C.mean():9.4f}  std = {C.std():8.4f}  n = {len(C)}')
    print()