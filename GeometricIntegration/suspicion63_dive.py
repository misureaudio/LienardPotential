"""
Suspicion 6.3, dive: where do the buggy and correct Newton IM orbits
diverge, and what does the buggy Newton do at that step?

Hypothesis: with the wrong-sign (2,1) Jacobian entry, Newton converges
(with residual below tolerance!) to a SPURIOUS root of the stage
equations (a cubic in u = xg - x), making the 'step' a wild jump that
is not the implicit midpoint step.

Run: D:/Source/hermes-dir/.venv/Scripts/python.exe suspicion63_dive.py
Writes: suspicion63_dive_results.txt
"""
import numpy as np

mu = 100.0
OUT = []
def p(s=""):
    OUT.append(s); print(s)

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

def stage_roots(x, v, h):
    """Real roots w of the stage equations written as a cubic in w = xg - x
    (endpoint increment; (mx,mv)=(x+w/2, w/h) is the midpoint):
      (h*mu/8) w^3 + (h*mu*x/2) w^2
      + (1 - (h*mu/2)(1-x^2) + h^2/4) w - h*v - (h^2*x/2) = 0
    (derived by substituting w = h*mv, vg = 2*mv - v into the two
    midpoint equations xg = x + h*mv, vg = v + h*f_v(mx,mv))."""
    c3 = h * mu / 8.0
    c2 = h * mu * x / 2.0
    c1 = 1.0 - (h * mu / 2.0) * (1.0 - x * x) + h * h / 4.0
    c0 = -h * v - h * h * x / 2.0
    r = np.roots([c3, c2, c1, c0])
    real = [float(z.real) for z in r if abs(z.imag) < 1e-9 * max(1.0, abs(z))]
    # residual check of each root in the ORIGINAL stage equations
    out = []
    for w in real:
        xg = x + w
        mv = w / h
        mx = (x + xg) / 2
        vg = 2 * mv - v
        res = max(abs(xg - x - h * mv), abs(vg - v - h * (mu * (1 - mx * mx) * mv - mx)))
        out.append((w, xg, vg, res))
    return out

for h in (0.5, 0.1, 0.01):
    p("=" * 74)
    p("h = %.2f" % h)
    p("=" * 74)
    xb, vb = 2.0, 0.0
    xc, vc = 2.0, 0.0
    found = False
    for i in range(int(3 * 162.8371 / h)):
        xgb, vgb, itb, resb, convb = buggy(xb, vb, h)
        xgc, vgc, itc, resc, convc = correct(xc, vc, h)
        dev = max(abs(xgb - xgc), abs(vgb - vgc))
        if dev > 1e-8:
            found = True
            p("  FIRST DIVERGENCE at step %d (t = %.3f)" % (i + 1, (i + 1) * h))
            p("  state before step: x = %.8f, v = %.8f" % (xb, vb))
            p("  correct step: xg = %.8f, vg = %.8f   (it=%d, res=%.1e, conv=%s)"
              % (xgc, vgc, itc, resc, convc))
            p("  buggy   step: xg = %.8f, vg = %.8f   (it=%d, res=%.1e, conv=%s)"
              % (xgb, vgb, itb, resb, convb))
            p("  deviation: %.3e" % dev)
            p("  ALL real roots of the stage equations at this state (u, xg, vg, residual):")
            for u, xg, vg, res in stage_roots(xb, vb, h):
                tag = ""
                if max(abs(u - (xgc - xb)), abs(vg - vgc)) < 1e-6:
                    tag = "  <-- correct Newton landed here (true IM step)"
                if max(abs(u - (xgb - xb)), abs(vg - vgb)) < 1e-6:
                    tag = "  <-- BUGGY Newton landed here (SPURIOUS root)"
                p("    u = %12.6f   xg = %12.6f   vg = %14.6f   res = %.1e%s"
                  % (u, xg, vg, res, tag))
            break
        xb, vb = xgb, vgb
        xc, vc = xgc, vgc
    if not found:
        p("  no divergence in 3 periods (orbits identical to 1e-8)")
    p("")

with open("suspicion63_dive_results.txt", "w") as fh:
    fh.write("\n".join(OUT) + "\n")
print("[written to suspicion63_dive_results.txt]")
