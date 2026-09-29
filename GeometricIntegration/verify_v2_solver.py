"""
Verify the CORRECTED robust implicit-midpoint solver reproduces the verified
ground truth (suspicion63_results.txt) BEFORE wiring it into the essay scripts.

Solver design (per the v2 remediation, Step 1):
  * exact Newton Jacobian  c = -h/2 * J21   (the sign the published scripts got wrong)
  * never silently return an unconverged iterate:
      - try Newton (tol=1e-12, itmax=200)
      - if it fails, fall back to an independent scipy.optimize.root (hybr)
      - take the smaller-residual candidate
      - if BOTH are garbage (residual > 1), raise RuntimeError (fail loud)
Run: ./.venv/Scripts/python.exe verify_v2_solver.py
"""
import numpy as np
from scipy.optimize import root

mu = 100.0
def fv(mx, mv):  return mu*(1.0-mx*mx)*mv - mx
def J21(mx, mv): return -2.0*mu*mx*mv - 1.0
def J22(mx, mv): return mu*(1.0-mx*mx)

def imid_step_robust(x, v, h, mu, tol=1e-12, itmax=200, fallback=True):
    """Robust implicit midpoint step (fail-loud). Returns (xg, vg)."""
    xg, vg = x, v
    res = None
    for _ in range(itmax):
        mx, mv = 0.5*(x+xg), 0.5*(v+vg)
        G1 = xg - x - h*mv
        G2 = vg - v - h*fv(mx, mv)
        r = max(abs(G1), abs(G2))
        res = (xg, vg, r)
        if r < tol:
            return xg, vg
        a, b = 1.0, -0.5*h
        c = -0.5*h*J21(mx, mv)          # EXACT (2,1) entry (published scripts had +)
        d = 1.0 - 0.5*h*J22(mx, mv)
        det = a*d - b*c
        xg -= ( d*G1 - b*G2)/det
        vg -= (-c*G1 + a*G2)/det
    if res is None or not (np.isfinite(res[0]) and np.isfinite(res[1])):
        raise RuntimeError("IM Newton diverged (non-finite) at x=%g v=%g h=%g" % (x, v, h))
    best = res
    if fallback:
        def F(z):
            mx, mv = 0.5*(x+z[0]), 0.5*(v+z[1])
            return [z[0]-x-h*mv, z[1]-v-h*(mu*(1-mx*mx)*mv-mx)]
        r = root(F, [x, v], method='hybr', tol=1e-12)
        cand = (float(r.x[0]), float(r.x[1]), max(abs(a) for a in F(r.x)))
        if np.isfinite(cand[0]) and np.isfinite(cand[1]) and cand[2] < best[2]:
            best = cand
    if best[2] > 1.0:
        raise RuntimeError("IM step UNRESOLVED (residual %.3e) at x=%g v=%g h=%g" % (best[2], x, v, h))
    return best[0], best[1]

def run(stepfn, h, t_end):
    x, v = 2.0, 0.0
    xprev = x; cross = []; mx, mv = abs(x), abs(v)
    for i in range(int(round(t_end/h))):
        x, v = stepfn(x, v, h)
        if not (np.isfinite(x) and np.isfinite(v)):
            return None, mx, mv, "DIVERGED"
        mx = max(mx, abs(x)); mv = max(mv, abs(v))
        if xprev < 0.0 and x >= 0.0 and v > 0.0: cross.append((i+1)*h)
        xprev = x
    T = float(np.mean(np.diff(cross)[-3:])) if len(cross) >= 2 else None
    return T, mx, mv, ("ok" if T else "no-crossings")

Tref = 162.8371
t3 = 3.0*Tref
print("Radau reference Tref = %.4f\n" % Tref)
print("  %-24s %10s %10s %10s  %s" % ("method (h)", "T", "max|x|", "max|v|", "status"))
for h in (0.5, 0.1, 0.01):
    T, mx, mv, st = run(lambda x,v,h: imid_step_robust(x,v,h,mu), h, t3)
    print("  %-24s %10s %10.3f %10.3f  %s" % ("IM robust, h=%.2f"%h,
          ("%.4f"%T) if T else "None", mx, mv, st))

print("\nGround truth (suspicion63_results.txt):")
print("  CORRECT IM h=0.50  T=None  max|x|=3.032  max|v|=15.94  (bounded, jump unresolved)")
print("  CORRECT IM h=0.10  T=286.5000  max|x|=2.954  max|v|=73.81")
print("  CORRECT IM h=0.01  T=179.0800  max|x|=2.090  max|v|=143.9")
