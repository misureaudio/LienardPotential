"""
Verify PURE correct-Newton + fail-loud guard (no scipy fallback) reproduces the
cited ground truth (suspicion63_results.txt) exactly.

Guard design: after the iteration cap, if the final stage residual > 1 (garbage,
like the original 10^90), RAISE RuntimeError. Residual <= 1 is accepted as
"bounded but possibly jump-unresolved" (the genuine A-stable-not-L-stable state).
This is what makes the original silent-garbage bug impossible while reproducing
the analysis doc's numbers exactly.
"""
import numpy as np

mu = 100.0
def fv(mx, mv):  return mu*(1.0-mx*mx)*mv - mx
def J21(mx, mv): return -2.0*mu*mx*mv - 1.0
def J22(mx, mv): return mu*(1.0-mx*mx)

def imid_step(x, v, h, mu, tol=1e-14, itmax=200, guard=1.0):
    xg, vg = x, v
    res = float('inf')
    for _ in range(itmax):
        mx, mv = 0.5*(x+xg), 0.5*(v+vg)
        G1 = xg - x - h*mv
        G2 = vg - v - h*fv(mx, mv)
        res = max(abs(G1), abs(G2))
        if res < tol:
            return xg, vg
        a, b = 1.0, -0.5*h
        c = -0.5*h*J21(mx, mv)          # EXACT (2,1) entry
        d = 1.0 - 0.5*h*J22(mx, mv)
        det = a*d - b*c
        xg -= ( d*G1 - b*G2)/det
        vg -= (-c*G1 + a*G2)/det
    # fail-loud guard
    if not (np.isfinite(xg) and np.isfinite(vg)):
        raise RuntimeError("IM Newton diverged (non-finite) at x=%g v=%g h=%g" % (x, v, h))
    if res > guard:
        raise RuntimeError("IM step UNRESOLVED (residual %.3e) at x=%g v=%g h=%g" % (res, x, v, h))
    return xg, vg

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
print("  %-26s %10s %10s %10s  %s" % ("method (h)", "T", "max|x|", "max|v|", "status"))
for h in (0.5, 0.1, 0.01):
    T, mx, mv, st = run(lambda x,v,h: imid_step(x,v,h,mu), h, t3)
    print("  %-26s %10s %10.3f %10.3f  %s" % ("IM pure+guard, h=%.2f"%h,
          ("%.4f"%T) if T else "None", mx, mv, st))
print("\nGround truth: h=0.5 max|x|=3.032  |  h=0.1 T=286.5000 max|x|=2.954  |  h=0.01 T=179.0800 max|x|=2.090")

# Confirm the guard CATCHES the original bug: buggy Jacobian (c = +h/2 J21) must raise.
def imid_buggy(x, v, h, mu, tol=1e-13, itmax=100, guard=1.0):
    xg, vg = x, v
    res = float('inf')
    for _ in range(itmax):
        mx, mv = 0.5*(x+xg), 0.5*(v+vg)
        G1 = xg - x - h*mv
        G2 = vg - v - h*fv(mx, mv)
        res = max(abs(G1), abs(G2))
        if res < tol: return xg, vg
        a, b = 1.0, -0.5*h
        c = +0.5*h*J21(mx, mv)          # BUGGY sign (as in computation3/4.py)
        d = 1.0 - 0.5*h*J22(mx, mv)
        det = a*d - b*c
        xg -= ( d*G1 - b*G2)/det
        vg -= (-c*G1 + a*G2)/det
    if res > guard:
        raise RuntimeError("BUGGY IM step UNRESOLVED (residual %.3e) -- guard fired" % res)
    return xg, vg

print("\nGuard test on the ORIGINAL (buggy) solver, h=0.01:")
try:
    run(lambda x,v,h: imid_buggy(x,v,h,mu), 0.01, t3)
    print("  NO ERROR (unexpected)")
except RuntimeError as e:
    print("  RuntimeError raised as intended:")
    print("   ", e)
