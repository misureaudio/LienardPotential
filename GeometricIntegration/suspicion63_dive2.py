"""
Suspicion 6.3, dive v2: track the growth of the deviation between the
buggy (published) and the correct Newton IM orbits, and characterize
the buggy Newton's iteration at the divergence point.
Also: verify the stage-equation cubic symbolically with sympy.

Run: D:/Source/hermes-dir/.venv/Scripts/python.exe suspicion63_dive2.py
Writes: suspicion63_dive2_results.txt
"""
import numpy as np
import sympy as sp

mu = 100.0
OUT = []
def p(s=""):
    OUT.append(s); print(s)

# ----------------------------------------------------------------------
p("=" * 74)
p("(A) Symbolic verification of the stage-equation cubic")
p("=" * 74)
x, v, h, mu_s = sp.symbols('x v h mu', real=True)
w = sp.symbols('w')   # w = xg - x, with mv = w/h, vg = 2*mv - v, mx = x + w/2
# stage equation 2: vg = v + h*( mu(1-mx^2) mv - mx )
lhs = 2 * w / h - v
rhs = v + h * (mu_s * (1 - (x + w / 2)**2) * (w / h) - (x + w / 2))
poly = sp.expand(lhs - rhs)
poly = sp.Poly(poly, w)
p("  stage equation 2 as polynomial in w (endpoint increment):")
for i in range(poly.degree(), -1, -1):
    p("    coeff w^%d = %s" % (i, sp.simplify(poly.coeff_monomial(w**i))))
# compare with the numeric formula (mv = w/h)
c3 = mu / 4.0
c2 = mu * x
c1 = h / 2.0 + mu * x * x - mu + 2.0 / h
c0 = h * x - 2.0 * v
sym3 = sp.simplify(poly.coeff_monomial(w**3) - c3)
sym2 = sp.simplify(poly.coeff_monomial(w**2) - c2)
sym1 = sp.simplify(poly.coeff_monomial(w**1) - c1)
sym0 = sp.simplify(poly.coeff_monomial(w**0) - c0)
p("  difference from numeric formula: c3 %s, c2 %s, c1 %s, c0 %s"
  % (sym3, sym2, sym1, sym0))
p("  (all zero => numeric cubic is correct)")
p("")

# ----------------------------------------------------------------------
def make_imid_step(buggy=False, tol=1e-14, itmax=200):
    sign = 1.0 if buggy else -1.0
    def step(x, v, h):
        xg, vg = x, v
        it = 0
        for it in range(itmax):
            mx, mv = 0.5 * (x + xg), 0.5 * (v + vg)
            G1 = xg - x - h * mv
            G2 = vg - v - h * (mu * (1 - mx * mx) * mv - mx)
            if abs(G1) < tol and abs(G2) < tol:
                return xg, vg, it, max(abs(G1), abs(G2)), True
            J21 = -2.0 * mu * mx * mv - 1.0
            J22 = mu * (1.0 - mx * mx)
            a, b = 1.0, -0.5 * h
            c = sign * 0.5 * h * J21
            d = 1.0 - 0.5 * h * J22
            det = a * d - b * c
            xg -= (d * G1 - b * G2) / det
            vg -= (-c * G1 + a * G2) / det
        return xg, vg, it, max(abs(G1), abs(G2)), False
    return step

correct = make_imid_step(buggy=False)
buggy = make_imid_step(buggy=True, tol=1e-13, itmax=100)

for h in (0.5, 0.1, 0.01):
    p("=" * 74)
    p("h = %.2f : deviation |buggy - correct| over time" % h)
    p("=" * 74)
    xb, vb = 2.0, 0.0
    xc, vc = 2.0, 0.0
    n = int(3 * 162.8371 / h)
    dev = 0.0
    marks = {1.0: None, 100.0: None, 1e6: None, 1e12: None}
    last_state = None
    for i in range(n):
        xgb, vgb, itb, resb, convb = buggy(xb, vb, h)
        xgc, vgc, itc, resc, convc = correct(xc, vc, h)
        d = max(abs(xgb - xgc), abs(vgb - vgc))
        if d > dev:
            dev = d
            for L, t in marks.items():
                if t is None and d > L:
                    marks[L] = (i + 1) * h
        if d > 1.0 and last_state is None:
            last_state = (i, xb, vb, xgb, vgb, itb, resb, convb, xgc, vgc, itc, resc, convc)
        xb, vb = xgb, vgb
        xc, vc = xgc, vgc
        if not (np.isfinite(xb) and np.isfinite(vb)):
            p("  buggy orbit left finite numbers at step %d" % (i + 1))
            break
    p("  |x|>1e6 first at: %s" % ("t=%.2f" % marks[1e6] if marks[1e6] else "-"))
    p("  max deviation over 3 periods: %.3e" % dev)
    if last_state:
        (i, xs, vs, xgb, vgb, itb, resb, convb, xgc, vgc, itc, resc, convc) = last_state
        p("  first step with deviation > 1: step %d (t = %.3f)" % (i + 1, (i + 1) * h))
        p("    state before:  x = %.8f, v = %.8f  (x-1 = %.4f)" % (xs, vs, xs - 1.0))
        p("    correct: xg = %.8f, vg = %.8f   it=%3d res=%.2e conv=%s"
          % (xgc, vgc, itc, resc, convc))
        p("    buggy  : xg = %.8g, vg = %.8g   it=%3d res=%.2e conv=%s"
          % (xgb, vgb, itb, resb, convb))
        # residual of the buggy iterate in the stage equations
        x, v = xs, vs
        xg, vg = xgb, vgb
        mx, mv = 0.5 * (x + xg), 0.5 * (v + vg)
        r1 = abs(xg - x - h * mv)
        r2 = abs(vg - v - h * (mu * (1 - mx * mx) * mv - mx))
        p("    buggy iterate stage residual (recomputed): G1 = %.3e, G2 = %.3e" % (r1, r2))
        # is the buggy step close to ANY root? solve cubic (mv = w/h)
        c3 = mu / 4.0
        c2 = mu * x
        c1 = h / 2.0 + mu * x * x - mu + 2.0 / h
        c0 = h * x - 2.0 * v
        rts = np.roots([c3, c2, c1, c0])
        real = [float(z.real) for z in rts if abs(z.imag) < 1e-8 * max(1.0, abs(z))]
        p("    real stage roots w at this state: %s"
          % ", ".join("%.4f" % w for w in real))
        if real:
            dmin = min(abs(xgb - (x + w)) for w in real)
            p("    |buggy xg - nearest root| = %.3e" % dmin)
        # where is the orbit: slow branch?
        if abs(xs) > 1.0:
            ys = xs / (mu * (1 - xs * xs))
            p("    (slow-manifold velocity at this x: y_s = %.4f; v = %.4f)" % (ys, vs))
    p("")

with open("suspicion63_dive2_results.txt", "w") as fh:
    fh.write("\n".join(OUT) + "\n")
print("[written to suspicion63_dive2_results.txt]")
