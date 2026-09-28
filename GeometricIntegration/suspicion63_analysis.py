"""
Suspicion 6.3 analysis (companion to geometric_integration_lienard_essay_v1.md).

Question: the essay's section 6.3 reports that the implicit midpoint rule
(symplectic + A-stable) on the mu=100 van der Pol cycle "escapes to ~1e30"
even at h=0.01. An independent re-implementation reports bounded orbits
(|x|,|v| up to ~100-150). Where do the ~1e30 numbers come from?

This script:
  (1) re-derives the EXACT Newton Jacobian for the implicit midpoint
      equations, symbolically (sympy), and compares it with the matrix
      actually used in geometric_integration_computation3/4.py;
  (2) implements the CORRECT implicit midpoint step (Newton, exact
      Jacobian) and cross-checks it against an independent solver
      (scipy.optimize.root);
  (3) runs the ORIGINAL (as-published) Newton and tracks, step by step,
      whether it actually converges to the root of the stage equations;
  (4) integrates the mu=100 van der Pol cycle with buggy vs correct
      implicit midpoint at h = 0.5, 0.1, 0.01, for 3 periods and longer;
  (5) reproduces the Strang-split and Radau rows of the essay's table
      (no Newton involved -> should reproduce exactly);
  (6) sanity checks on the harmonic oscillator.

Run: D:/Source/hermes-dir/.venv/Scripts/python.exe suspicion63_analysis.py
Writes: suspicion63_results.txt, suspicion63_figure.png
"""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import root

OUT = []
def p(s=""):
    OUT.append(s)
    print(s)

mu = 100.0

# ----------------------------------------------------------------------
p("=" * 74)
p("(1) SYMBOLIC: exact Newton Jacobian for the implicit midpoint stage")
p("=" * 74)
x, v, h, mu_s = sp.symbols('x v h mu', real=True)
xg, vg = sp.symbols('xg vg', real=True)
mx = (x + xg) / 2
mv = (v + vg) / 2
fx = mv
fv = mu_s * (1 - mx**2) * mv - mx
G1 = xg - x - h * fx
G2 = vg - v - h * fv
JG = sp.Matrix([[sp.diff(G1, xg), sp.diff(G1, vg)],
                [sp.diff(G2, xg), sp.diff(G2, vg)]])
p("Exact Jacobian dG/d(xg,vg) of the stage equations G = 0:")
p(sp.pretty(JG))
p("  entry (2,1) = dG2/dxg = -h/2 * J21,  with J21 = dfv/dmx = -2*mu*mx*mv - 1")
J21 = -2*mu_s*mx*mv - 1
p("  i.e.  c_exact = %s" % sp.simplify(-sp.Rational(1,2)*h*J21))
p("The scripts computation3/4.py use  c_script = +h/2*J21 = %s" % sp.simplify(sp.Rational(1,2)*h*J21))
p("=> the (2,1) entry has the WRONG SIGN in the published scripts.")
p("   (computation2.py had c = -h/2*J21 [correct] but the v-update line")
p("    vg -= (c*G1 - a*G2)/det [wrong sign]; the 'fix' in 3/4 flipped c")
p("    instead of the update line, so the effective Jacobian is still wrong.)")
p("")
# Newton map for harmonic oscillator with the CORRECT Jacobian
xs, vs = sp.symbols('xs vs')
eq1 = sp.Eq(xs, x + h*(v + vs)/2)
eq2 = sp.Eq(vs, v - h*(x + xs)/2)
sol = sp.solve([eq1, eq2], [xs, vs], dict=True)[0]
xs_s = sp.simplify(sol[xs]); vs_s = sp.simplify(sol[vs])
M = sp.Matrix([[sp.diff(xs_s, x), sp.diff(xs_s, v)],
               [sp.diff(vs_s, x), sp.diff(vs_s, v)]])
J = sp.Matrix([[0, 1], [-1, 0]])
p("Harmonic oscillator: exact IM endpoint map (closed form)")
p("  x1 = %s" % xs_s)
p("  v1 = %s" % vs_s)
p("  det M = %s ;  M^T J M - J = %s   (symplectic, as essay eq. (8))"
  % (sp.simplify(M.det()), sp.simplify(M.T * J * M - J)))
p("")

# ----------------------------------------------------------------------
p("=" * 74)
p("(2) NUMERICAL: correct IM step vs independent solver (scipy.optimize.root)")
p("=" * 74)

def vdp_f(x, v, mu):
    return (v, mu * (1 - x * x) * v - x)

def make_imid_step(fv, J21fn, J22fn, buggy=False, tol=1e-14, itmax=200):
    """Implicit midpoint step for f = (v, fv(mx,mv)), Newton on the stage
    equations. buggy=True reproduces the published scripts
    (c = +h/2*J21 instead of the exact c = -h/2*J21)."""
    sign = 1.0 if buggy else -1.0
    def step(x, v, h):
        xg, vg = x, v
        it = 0
        for it in range(itmax):
            mx, mv = 0.5 * (x + xg), 0.5 * (v + vg)
            G1 = xg - x - h * mv
            G2 = vg - v - h * fv(mx, mv)
            if abs(G1) < tol and abs(G2) < tol:
                return xg, vg, it, max(abs(G1), abs(G2)), True
            a, b = 1.0, -0.5 * h
            c = sign * 0.5 * h * J21fn(mx, mv)
            d = 1.0 - 0.5 * h * J22fn(mx, mv)
            det = a * d - b * c
            xg -= (d * G1 - b * G2) / det
            vg -= (-c * G1 + a * G2) / det
        return xg, vg, it, max(abs(G1), abs(G2)), False
    return step

vdp_fv = lambda mx, mv: mu * (1.0 - mx * mx) * mv - mx
vdp_J21 = lambda mx, mv: -2.0 * mu * mx * mv - 1.0
vdp_J22 = lambda mx, mv: mu * (1.0 - mx * mx)

# correct Newton (exact Jacobian) vs the Newton as published (computation3/4.py)
imid_step_correct = make_imid_step(vdp_fv, vdp_J21, vdp_J22, buggy=False)
imid_step_buggy = make_imid_step(vdp_fv, vdp_J21, vdp_J22, buggy=True,
                                 tol=1e-13, itmax=100)

def imid_root_indep(x, v, h, mu):
    """Independent solve of the stage equations (scipy hybr, no hand Jacobian)."""
    def F(z):
        mx, mv = 0.5 * (x + z[0]), 0.5 * (v + z[1])
        return [z[0] - x - h * mv,
                z[1] - v - h * (mu * (1 - mx * mx) * mv - mx)]
    r = root(F, [x, v], method='hybr', tol=1e-14)
    return r.x, r.success, max(abs(a) for a in F(r.x))

# cross-check at a few states (mild and stiff, but within Newton's basin)
states = [(2.0, 0.0, 0.01), (1.5, -30.0, 0.01), (0.3, -95.0, 0.01),
          (2.0, 0.0, 0.5), (1.5, -20.0, 0.5)]
p("  state (x,v)          h     Newton(correct)        scipy.root         agree")
ok = True
for (x0, v0, h0) in states:
    xn, vn, it, res, conv = imid_step_correct(x0, v0, h0)
    (xr, vr), s2, resr = imid_root_indep(x0, v0, h0, mu)
    err = max(abs(xn - xr), abs(vn - vr))
    ok = ok and conv and s2 and err < 1e-8
    p("  (%6.2f,%8.2f)   %5.2f   it=%2d res=%.1e      res=%.1e        %.2e %s"
      % (x0, v0, h0, it, res, resr, err, "OK" if err < 1e-8 else "MISMATCH"))
p("  => correct Newton == independent root solve: %s" % ("YES" if ok else "NO"))
p("")

# ----------------------------------------------------------------------
p("=" * 74)
p("(3) Does the PUBLISHED Newton actually converge? (mu=100)")
p("    Per-step diagnostics: iterations used, final residual, converged?")
p("=" * 74)
for h0 in (0.5, 0.1, 0.01):
    for name, stpf in (("correct", imid_step_correct), ("buggy  ", imid_step_buggy)):
        x, v = 2.0, 0.0
        nconv = 0
        iters = []
        fails = []
        for i in range(120):
            x, v, it, res, conv = stpf(x, v, h0)
            iters.append(it)
            if conv:
                nconv += 1
            else:
                fails.append((i + 1, it, res))
            if not (np.isfinite(x) and np.isfinite(v)):
                break
        p("  h=%5.2f %s: converged %3d/120 steps;  max iters %3d"
          % (h0, name, nconv, max(iters)))
        if fails:
            f = fails[0]
            p("             first non-converged step: %d (iters=%d, residual=%.2e);  total failed: %d"
              % (f[0], f[1], f[2], len(fails)))
p("")
p("  A Newton step that converges lands on a root of G=0 -- but with the")
p("  WRONG-sign Jacobian it can lock onto a spurious root of the (cubic in")
p("  the stage velocity) stage system, far from the nearby true-IM root.")
p("  Then the 'step' is not the implicit midpoint step at all, even though")
p("  it satisfies the stage equations.")
p("")
p("  Step-by-step comparison over 120 steps from (2,0):")
for h0 in (0.5, 0.1, 0.01):
    xb, vb, xc, vc = 2.0, 0.0, 2.0, 0.0
    maxdev, bigjumps = 0.0, []
    for i in range(120):
        xcb, vcb = imid_step_buggy(xb, vb, h0)[:2]
        xcc, vcc = imid_step_correct(xc, vc, h0)[:2]
        dev = max(abs(xcb - xcc), abs(vcb - vcc))
        maxdev = max(maxdev, dev)
        if dev > 1.0:
            bigjumps.append((i + 1, dev))
        xb, vb, xc, vc = xcb, vcb, xcc, vcc
    p("  h=%5.2f: max |buggy - correct| step = %.3e ;  jumps >1 at steps: %s"
      % (h0, maxdev, (bigjumps[:5] + ["..."]) if len(bigjumps) > 5 else bigjumps))
p("")

# ----------------------------------------------------------------------
p("=" * 74)
p("(4) Orbits on the mu=100 van der Pol cycle, start (2,0)")
p("    True: T=162.837, amplitude 2, max|v| ~ 134")
p("=" * 74)

def ref_period(mu, rtol=1e-10):
    def f(t, y):
        x, v = y
        return [v, mu * (1 - x * x) * v - x]
    def cross(t, y):
        return y[0]
    cross.terminal = False
    cross.direction = 1
    sol = solve_ivp(f, (0.0, 4.0 * mu * 1.614 + 5.0), [2.0, 0.0], method="Radau",
                    rtol=rtol, atol=rtol, events=cross, max_step=0.2)
    per = np.diff(sol.t_events[0])
    return float(np.mean(per[-3:]))

Tref = ref_period(mu)
p("  Radau reference period Tref = %.4f" % Tref)
p("")

def run_orbit(stepfn, h, mu, t_end, track=False):
    x, v = 2.0, 0.0
    xprev = x
    cross = []
    mx, mv = abs(x), abs(v)
    first_big = {L: None for L in (100.0, 1e6, 1e12, 1e20)}
    n = int(round(t_end / h))
    for i in range(n):
        x, v = stepfn(x, v, h)
        if not (np.isfinite(x) and np.isfinite(v)):
            return dict(T=None, mx=mx, mv=mv, n=i + 1, status="DIVERGED",
                        cross=len(cross), first_big=first_big)
        if abs(x) > mx: mx = abs(x)
        if abs(v) > mv: mv = abs(v)
        for L in first_big:
            if first_big[L] is None and abs(x) > L:
                first_big[L] = (i + 1) * h
        if xprev < 0.0 and x >= 0.0 and v > 0.0:
            cross.append((i + 1) * h)
        xprev = x
    T = float(np.mean(np.diff(cross)[-3:])) if len(cross) >= 2 else None
    return dict(T=T, mx=mx, mv=mv, n=n,
                status=("ok" if T else "no-crossings"), cross=len(cross),
                first_big=first_big, track=track)

def imid_b(x, v, h):
    return imid_step_buggy(x, v, h)[:2]

def imid_c(x, v, h):
    return imid_step_correct(x, v, h)[:2]

t3 = 3.0 * Tref
p("  --- 3 periods (~%.1f time units) ---" % t3)
p("  %-28s %12s %12s %10s %12s  %s" % ("method (h)", "T", "max|x|", "max|v|", "crossings", "|x|>1e6 at t="))
for h0 in (0.5, 0.1, 0.01):
    r = run_orbit(imid_b, h0, mu, t3)
    tbig = r["first_big"][1e6]
    p("  %-28s %12s %12.3e %10.3e %9d   %s"
      % ("BUGGY IM (published), h=%.2f" % h0,
         ("%.4f" % r["T"]) if r["T"] else "None", r["mx"], r["mv"], r["cross"],
         ("t=%.2f" % tbig) if tbig else "-"))
for h0 in (0.5, 0.1, 0.01):
    r = run_orbit(imid_c, h0, mu, t3)
    tbig = r["first_big"][1e6]
    p("  %-28s %12s %12.3e %10.3e %9d   %s"
      % ("CORRECT IM, h=%.2f" % h0,
         ("%.4f" % r["T"]) if r["T"] else "None", r["mx"], r["mv"], r["cross"],
         ("t=%.2f" % tbig) if tbig else "-"))
p("")

# longer runs, correct method: does it stay bounded?
p("  --- longer runs, CORRECT implicit midpoint ---")
p("  %-26s %12s %12s %10s %10s" % ("h, periods", "T (last 3)", "max|x|", "max|v|", "status"))
longer = [(0.01, 20), (0.1, 10), (0.5, 5)]
for h0, nper in longer:
    r = run_orbit(imid_c, h0, mu, nper * Tref)
    p("  %-26s %12s %12.3e %10.3e %10s"
      % ("h=%.2f, %d per" % (h0, nper),
         ("%.4f" % r["T"]) if r["T"] else "None", r["mx"], r["mv"], r["status"]))
p("")

# per-period spike (max|v| within each period) for the correct method at h=0.01
p("  --- per-period max|v| spikes, CORRECT IM, h=0.01, 10 periods (ringing check) ---")
x, v = 2.0, 0.0
h0 = 0.01
per_steps = int(round(Tref / h0))
spikes = []
for k in range(10):
    mx, mv = 0.0, 0.0
    for i in range(per_steps):
        x, v = imid_c(x, v, h0)
        mx = max(mx, abs(x)); mv = max(mv, abs(v))
    spikes.append((k + 1, mx, mv))
for k, mx, mv in spikes:
    p("    period %2d: max|x| = %9.3f   max|v| = %9.3f" % (k, mx, mv))
p("")

# ----------------------------------------------------------------------
p("=" * 74)
p("(5) Strang split (exact damping) -- NO Newton involved, must reproduce")
p("=" * 74)

def strang_split(x, v, h, mu):
    v *= np.exp(mu * (1 - x * x) * h / 2.0)
    v1 = v - 0.5 * h * x
    x1 = x + h * v1
    v2 = v1 - 0.5 * h * x1
    v2 *= np.exp(mu * (1 - x1 * x1) * h / 2.0)
    return x1, v2

for h0 in (0.5, 0.2):
    r = run_orbit(lambda x, v, h: strang_split(x, v, h, mu), h0, mu, t3)
    p("  Strang split, h=%.2f : T=%s  max|x|=%.3e  max|v|=%.3e  %s"
      % (h0, ("%.4f" % r["T"]) if r["T"] else "None", r["mx"], r["mv"], r["status"]))
p("  (essay table: h=0.5 -> ~2.0e15 blow-up; h=0.2 -> T=151.40, max|x|~4.6e5)")
p("")

# ----------------------------------------------------------------------
p("=" * 74)
p("(6) Radau IIA rows of the essay table (L-stable, not symplectic)")
p("=" * 74)
def radau_period(h, mu, t_end):
    def f(t, y):
        x, v = y
        return [v, mu * (1 - x * x) * v - x]
    def cross(t, y):
        return y[0]
    cross.terminal = False
    cross.direction = 1
    sol = solve_ivp(f, (0.0, t_end), [2.0, 0.0], method="Radau", rtol=1e-9, atol=1e-9,
                    events=cross, max_step=h)
    if len(sol.t_events[0]) < 2:
        return None
    return float(np.mean(np.diff(sol.t_events[0])[-3:]))
for h0 in (0.5, 1.0):
    T = radau_period(h0, mu, 3 * Tref)
    p("  Radau IIA max_step=%.1f : T=%s  (ref %.4f)" % (h0, ("%.4f" % T) if T else "None", Tref))
p("")

# ----------------------------------------------------------------------
p("=" * 74)
p("(7) Stability-function values at the relevant z = lambda*h")
p("=" * 74)
def R(z):
    return (1 + z / 2) / (1 - z / 2)
rows = [
    ("fast mode at x=2, lambda=-300", -300 * 0.01, "h=0.01"),
    ("fast mode at x=2, lambda=-300", -300 * 0.1, "h=0.1"),
    ("fast mode at x=2, lambda=-300", -300 * 0.5, "h=0.5"),
    ("central region x=0, lambda=+100", 100 * 0.01, "h=0.01"),
    ("central region x=0, lambda=+100", 100 * 0.1, "h=0.1"),
    ("central region x=0, lambda=+100", 100 * 0.5, "h=0.5"),
]
p("  %-34s %10s %8s | R(z) %10s | e^z %12s" % ("mode", "z", "h", "", ""))
for label, z, tag in rows:
    p("  %-34s %10.2f %8s | %10.5f | %12.3e" % (label, z, tag, R(z), np.exp(z) if z < 700 else float('inf')))
p("  R(-inf) = -1  (A-stable, not L-stable): fast modes are damped to |R| ~ 1,")
p("  with a sign flip each step, instead of to 0.")
p("")

# ----------------------------------------------------------------------
p("=" * 74)
p("(8) Sanity: harmonic oscillator, 1e5 steps, h=0.1 (mild -> Newton robust)")
p("=" * 74)
harm_fv = lambda mx, mv: -mx
harm_J21 = lambda mx, mv: -1.0
harm_J22 = lambda mx, mv: 0.0
imid_c_h = make_imid_step(harm_fv, harm_J21, harm_J22, buggy=False)
imid_b_h = make_imid_step(harm_fv, harm_J21, harm_J22, buggy=True,
                          tol=1e-13, itmax=100)
for name, fn in (("correct", imid_c_h), ("buggy  ", imid_b_h)):
    x, v = 1.0, 0.0
    E0 = 0.5
    mxe = 0.0
    for i in range(100000):
        x, v, _, _, _ = fn(x, v, 0.1)
        mxe = max(mxe, abs(0.5 * x * x + 0.5 * v * v - E0))
    p("  %s IM: max|E-E0| = %.3e  (bounded: no drift)" % (name, mxe))
p("")

with open("suspicion63_results.txt", "w") as fh:
    fh.write("\n".join(OUT) + "\n")
p("[written to suspicion63_results.txt]")

# ----------------------------------------------------------------------
# Figure: buggy vs correct IM, h=0.01, mu=100, ~3.5 periods
# ----------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def collect(fn, h, t_end):
    x, v = 2.0, 0.0
    ts, xs, vs = [0.0], [2.0], [0.0]
    n = int(round(t_end / h))
    for i in range(n):
        x, v = fn(x, v, h)
        ts.append((i + 1) * h); xs.append(x); vs.append(v)
    return np.array(ts), np.array(xs), np.array(vs)

t_end = 3.5 * Tref
tb, xb, vb = collect(imid_b, 0.01, t_end)
tc, xc, vc = collect(imid_c, 0.01, t_end)

fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
ax[0].plot(tb, np.abs(xb), lw=0.6, color="tab:red")
ax[0].set_yscale("log")
ax[0].set_ylabel("|x| (log)")
ax[0].set_title("mu=100 van der Pol, h=0.01: implicit midpoint\n"
                "RED: published Newton (wrong Jacobian sign)   BLUE: correct Newton")
ax[0].grid(alpha=0.3)
ax[1].plot(tc, xc, lw=0.7, color="tab:blue", label="correct IM: x(t)")
ax[1].set_ylabel("x(t)")
ax[1].set_xlabel("t")
ax[1].legend()
ax[1].grid(alpha=0.3)
fig.tight_layout()
fig.savefig("suspicion63_figure.png", dpi=140)
p("[figure saved: suspicion63_figure.png]")
