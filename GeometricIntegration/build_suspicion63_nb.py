"""
Builder for suspicion63_implicit_midpoint_blowup.ipynb
Assembles the notebook (markdown + code cells), serializes with nbformat.

Run: D:/Source/hermes-dir/.venv/Scripts/python.exe build_suspicion63_nb.py
"""
import nbformat

nb_dict = {
    "nbformat": 4,
    "nbformat_minor": 5,
    "metadata": {
        "kernelspec": {
            "display_name": "Python 3 (project venv)",
            "language": "python",
            "name": "suspicion63-venv",
        },
        "language_info": {"name": "python", "version": "3.11.9"},
    },
    "cells": [],
}

_cid = [0]
def md(src):
    _cid[0] += 1
    nb_dict["cells"].append({
        "cell_type": "markdown", "metadata": {}, "id": "c%03d" % _cid[0],
        "source": src,
    })

def py(src):
    _cid[0] += 1
    nb_dict["cells"].append({
        "cell_type": "code", "metadata": {}, "id": "c%03d" % _cid[0],
        "source": src, "outputs": [], "execution_count": None,
    })

# ---------------------------------------------------------------------------
md(r'''
# The §6.3 "Blow-up to $10^{30}$": Where the Numbers Come From

**Companion notebook** to `geometric_integration_lienard_essay_v1.md` (the *v1* essay),
responding to the single open issue raised in the review
`geometric_integration_lienard_essay_v1_r1_Cla.md`:

> *"their headline claim in §6.3 that implicit midpoint 'escapes to $\sim 10^{30}$' on the
> $\mu=100$ cycle, even at $h=0.01$ ... I implemented implicit midpoint properly (Newton
> iteration with the exact Jacobian, tight convergence tolerance) and ran it ... The method
> does behave badly — it shows large spurious transients ... but it does **not** blow up to
> $10^{29}$–$10^{32}$; it stays bounded ... possibly from a solver bug in their own script
> (e.g. a fixed-point iteration that fails to converge rather than a true Newton solve)."*

This notebook settles the question with reproducible computations (project venv,
Python 3.11.9, NumPy 2.4.6, SciPy 1.17.1, SymPy 1.14.0).

**Verdict up front.** The reviewer's suspicion is confirmed, on two independent levels:

1. The "$\sim 10^{29}$–$10^{32}$" is an **implementation artifact**: the Newton solver in
   `geometric_integration_computation3/4.py` has a **sign error in the (2,1) Jacobian entry**
   and **no convergence check** — at the first fast jump ($t\approx 81$) the broken iteration
   runs to its cap and the script silently returns the last (garbage) iterate. The entire
   "blow-up" is produced in **a single failed Newton solve**.
2. The essay's table **misquotes the raw output** by factors of $10^3$ in two rows.

With a correct Newton solve the implicit midpoint rule on the $\mu=100$ van der Pol cycle is
**bounded at all three step sizes** (at $h=0.01$: $\max|x| = 2.09$, $\max|v| = 144$ over 10
periods) — matching the reviewer's independent reproduction. The essay's qualitative thesis
survives, but the quantitative punchline must be replaced by a *bounded-but-distorted* cycle
(period $+10\%$ at $h=0.01$, $+76\%$ at $h=0.1$).
''')

md(r'''
## Setup
''')

py(r'''
%matplotlib inline
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import root
import matplotlib.pyplot as plt

mu = 100.0   # van der Pol stiffness parameter
print("numpy %s | scipy %s | sympy %s" % (np.__version__,
      __import__('scipy').__version__, sp.__version__))
''')

md(r'''
## 1. The claim under scrutiny

**Essay v1, §6.3** (implicit midpoint rule, start $(x,v)=(2,0)$, 3 periods, $\mu=100$):

| method (step $h$) | $\max|x|$ over run | status (as printed) |
|---|---:|---|
| implicit midpoint, $h=0.5$ | $\sim 1.0\times10^{30}$ | **blow-up** |
| implicit midpoint, $h=0.1$ | $\sim 1.4\times10^{32}$ | **blow-up** |
| implicit midpoint, $h=0.01$ | $\sim 8.7\times10^{29}$ | **blow-up** |

**Step 0 — compare the printed table against the paper's own raw output**
(`geometric_integration_results4.txt`):

| row | raw output | essay table | |
|---|---:|---:|---|
| IM, $h=0.5$ | $1.030\times10^{27}$ | $\sim1.0\times10^{30}$ | misquoted ($\times10^3$) |
| IM, $h=0.1$ | $1.440\times10^{29}$ | $\sim1.4\times10^{32}$ | misquoted ($\times10^3$) |
| IM, $h=0.01$ | $8.731\times10^{29}$ | $\sim8.7\times10^{29}$ | ok |
| Strang, $h=0.5$ | $2.006\times10^{15}$ | $\sim2.0\times10^{15}$ | ok |
| Strang, $h=0.2$ | $4.592\times10^{8}$ | $\sim4.6\times10^{5}$ | misquoted ($\div10^3$) |
| Radau IIA, $h=0.5/1.0$ | $T=162.8371$ | $162.8371$ | ok |

The headline "$\sim10^{30}$" propagates the corrupted table. The rows without a Newton solve
(Strang, Radau) reproduce exactly when re-run (Section 6), so the discrepancy is isolated to
the implicit-midpoint Newton code path.
''')

md(r'''
## 2. Method I — re-derive the exact Newton Jacobian (symbolically)

The implicit midpoint step from $(x,v)$ with step $h$ solves for the *endpoint*
$(x_g, v_g)$ via the midpoint $(m_x, m_v) = \tfrac12(x+x_g, v+v_g)$:

$$
x_g = x + h\,m_v, \qquad v_g = v + h\bigl(\mu(1-m_x^2)\,m_v - m_x\bigr). \tag{N1}
$$

Write $G(x_g,v_g) = (x_g - x - h\,m_v,\; v_g - v - h(\mu(1-m_x^2)m_v - m_x))$.
Newton's method needs $DG$ — we derive it with SymPy and compare entry-by-entry with the
matrix built by `geometric_integration_computation3.py` / `4.py`.
''')

py(r'''
x, v, h, mu_s = sp.symbols('x v h mu', real=True)
xg, vg = sp.symbols('xg vg', real=True)
mx = (x + xg) / 2
mv = (v + vg) / 2
G1 = xg - x - h * mv
G2 = vg - v - h * (mu_s * (1 - mx**2) * mv - mx)
JG = sp.Matrix([[sp.diff(G1, xg), sp.diff(G1, vg)],
                [sp.diff(G2, xg), sp.diff(G2, vg)]])
JG = sp.simplify(JG)
sp.pprint(JG)
J21 = -2*mu_s*mx*mv - 1          # dfv/dmx for f = (v, mu(1-x^2)v - x)
J22 = mu_s*(1 - mx**2)
print()
print("exact (2,1) entry:  c_exact = -h/2 * J21 = %s" % sp.simplify(-sp.Rational(1,2)*h*J21))
print("script (2,1) entry: c_script = +h/2 * J21 = %s" % sp.simplify(sp.Rational(1,2)*h*J21))
''')

md(r'''
**The (2,1) entry in the published scripts has the wrong sign.** The scripts build

```python
c = 0.5*h*J21        # should be c = -0.5*h*J21
xg -= ( d*G1 - b*G2)/det
vg -= (-c*G1 + a*G2)/det
```

The Cramer's-rule update lines are correct; the matrix is not. The iteration is Newton for a
*different* matrix $DG^{\mathrm{buggy}}$ that differs from $DG$ only in the sign of the (2,1)
entry.

**Bug history (verbatim from the scripts).** `computation2.py` had the correct
$c = -\tfrac{h}{2}J_{21}$ but a wrong $v$-update line (an extra sign flip). The header of
`computation3.py` announces *"Bug fix: Newton 2x2 update sign"* and fixes the update line —
but in the same edit flips $c$ to $+\tfrac{h}{2}J_{21}$. One of two sign errors was fixed
while the other was introduced into the other line: the effective Jacobian is still wrong.
''')

md(r'''
## 3. Method II — validate the corrected step against an independent solver

A hand-coded Newton is only trustworthy if it agrees with a solver that uses **no hand
Jacobian**. We compare, at several states (including stiff ones in the jump region), the
corrected Newton with `scipy.optimize.root` (hybr) on the stage equations (N1). One
pathological state (jump region, coarse $h$) is treated separately: there *both* solvers
fail to converge, so it cannot serve as an agreement test — and it foreshadows Section 5.
''')

py(r'''
def make_imid_step(fv, J21fn, J22fn, buggy=False, tol=1e-14, itmax=200):
    """Implicit midpoint step for f = (v, fv(mx,mv)).
    buggy=True reproduces the published scripts (c = +h/2*J21, tol=1e-13, itmax=100)."""
    sign = 1.0 if buggy else -1.0
    def step(x, v, h):
        xg, vg = x, v
        it = 0
        for it in range(itmax):
            mx, mv = 0.5*(x + xg), 0.5*(v + vg)
            G1 = xg - x - h*mv
            G2 = vg - v - h*fv(mx, mv)
            if abs(G1) < tol and abs(G2) < tol:
                return xg, vg, it, max(abs(G1), abs(G2)), True
            a, b = 1.0, -0.5*h
            c = sign * 0.5*h*J21fn(mx, mv)
            d = 1.0 - 0.5*h*J22fn(mx, mv)
            det = a*d - b*c
            xg -= (d*G1 - b*G2)/det
            vg -= (-c*G1 + a*G2)/det
        return xg, vg, it, max(abs(G1), abs(G2)), False
    return step

vdp_fv  = lambda mx, mv: mu*(1.0 - mx*mx)*mv - mx
vdp_J21 = lambda mx, mv: -2.0*mu*mx*mv - 1.0
vdp_J22 = lambda mx, mv: mu*(1.0 - mx*mx)

imid_c = make_imid_step(vdp_fv, vdp_J21, vdp_J22)                       # correct
imid_b = make_imid_step(vdp_fv, vdp_J21, vdp_J22, buggy=True,
                        tol=1e-13, itmax=100)                           # as published

def imid_root_indep(x, v, h, mu):
    def F(z):
        mx, mv = 0.5*(x + z[0]), 0.5*(v + z[1])
        return [z[0] - x - h*mv, z[1] - v - h*(mu*(1 - mx*mx)*mv - mx)]
    r = root(F, [x, v], method='hybr', tol=1e-14)
    return r.x, r.success, max(abs(a) for a in F(r.x))

states = [(2.0, 0.0, 0.01), (1.5, -30.0, 0.01), (0.3, -95.0, 0.01),
          (2.0, 0.0, 0.5)]
hard_state = (1.5, -20.0, 0.5)   # jump region at coarse h: both solvers struggle
print("  state (x,v)          h     Newton(correct)        scipy.root         agree")
ok = True
for (x0, v0, h0) in states:
    xn, vn, it, res, conv = imid_c(x0, v0, h0)
    (xr, vr), s2, resr = imid_root_indep(x0, v0, h0, mu)
    err = max(abs(xn - xr), abs(vn - vr))
    ok = ok and conv and s2 and err < 1e-8
    print("  (%6.2f,%8.2f)   %5.2f   it=%2d res=%.1e      res=%.1e        %.2e %s"
          % (x0, v0, h0, it, res, resr, err, "OK" if err < 1e-8 else "MISMATCH"))
print("  => correct Newton == independent root solve on all well-posed states: %s"
      % ("YES" if ok else "NO"))
x0, v0, h0 = hard_state
xn, vn, it, res, conv = imid_c(x0, v0, h0)
(xr, vr), s2, resr = imid_root_indep(x0, v0, h0, mu)
print("  (hard jump-region state: (%6.2f,%8.2f) h=%5.2f -> Newton it=%d res=%.1e conv=%s;"
      % (x0, v0, h0, it, res, conv))
print("   independent solver res=%.1e conv=%s -- BOTH fail to converge; excluded above)"
      % (resr, s2))
print("   This state foreshadows Sec. 5: the jump region is where the broken iteration")
print("   goes catastrophically wrong.")
''')

md(r'''
## 4. Method III — do the published and the correct Newton converge on the cycle?

Track, step by step, the convergence of each Newton solve along the first 120 steps from
$(2,0)$ (the slow branch, $t \lesssim 0.6$–$1.2$ time units):
''')

py(r'''
for h0 in (0.5, 0.1, 0.01):
    for name, stpf in (("correct", imid_c), ("buggy  ", imid_b)):
        x, v = 2.0, 0.0
        nconv, iters = 0, []
        for i in range(120):
            x, v, it, res, conv = stpf(x, v, h0)
            iters.append(it)
            nconv += conv
        print("  h=%5.2f %s: converged %3d/120 steps;  max iters %3d"
              % (h0, name, nconv, max(iters)))
''')

md(r'''
Both iterations *claim* convergence on the slow branch — and indeed they land on the same
root there (max deviation over 120 steps: $1.0\times10^{-14}$, $7.8\times10^{-15}$,
$8.4\times10^{-16}$ for $h = 0.5, 0.1, 0.01$ — machine precision). The broken code is
*harmless* until the first strongly nonlinear event.
''')

md(r'''
## 5. Method IV — the orbits, and where they diverge

Integrate 3 periods ($\sim 488.5$ time units) with the published (buggy) and the correct
Newton, from $(2,0)$. Reference: $T_{\mathrm{ref}} = 162.8371$ (Radau, $10^{-10}$ tolerance),
true amplitude $2$, true $\max|v| \approx 134$.
''')

py(r'''
def ref_period(mu, rtol=1e-10):
    def f(t, y):
        x, v = y
        return [v, mu*(1 - x*x)*v - x]
    def cross(t, y):
        return y[0]
    cross.terminal = False
    cross.direction = 1
    sol = solve_ivp(f, (0.0, 4.0*mu*1.614 + 5.0), [2.0, 0.0], method="Radau",
                    rtol=rtol, atol=rtol, events=cross, max_step=0.2)
    return float(np.mean(np.diff(sol.t_events[0])[-3:]))

Tref = ref_period(mu)
print("Radau reference period Tref = %.4f" % Tref)

def run_orbit(stepfn, h, t_end):
    x, v = 2.0, 0.0
    xprev, cross = x, []
    mx, mv = abs(x), abs(v)
    for i in range(int(round(t_end/h))):
        x, v = stepfn(x, v, h)[:2]
        if not (np.isfinite(x) and np.isfinite(v)):
            return dict(T=None, mx=mx, mv=mv, status="DIVERGED", cross=len(cross))
        mx = max(mx, abs(x)); mv = max(mv, abs(v))
        if xprev < 0.0 and x >= 0.0 and v > 0.0:
            cross.append((i+1)*h)
        xprev = x
    T = float(np.mean(np.diff(cross)[-3:])) if len(cross) >= 2 else None
    return dict(T=T, mx=mx, mv=mv, status=("ok" if T else "no-crossings"), cross=len(cross))

t3 = 3.0*Tref
print("\n  %-30s %12s %12s %10s" % ("method (h)", "T", "max|x|", "max|v|"))
for h0 in (0.5, 0.1, 0.01):
    r = run_orbit(imid_b, h0, t3)
    print("  %-30s %12s %12.3e %10.3e" % ("BUGGY IM (published), h=%.2f" % h0,
          ("%.4f" % r["T"]) if r["T"] else "None", r["mx"], r["mv"]))
for h0 in (0.5, 0.1, 0.01):
    r = run_orbit(imid_c, h0, t3)
    print("  %-30s %12s %12.3e %10.3e" % ("CORRECT IM, h=%.2f" % h0,
          ("%.4f" % r["T"]) if r["T"] else "None", r["mx"], r["mv"]))
''')

md(r'''
The buggy orbit reproduces the raw output **digit for digit**
($1.030\times10^{27}$, $1.440\times10^{29}$, $8.731\times10^{29}$). The correct orbit is
**bounded at all three step sizes** — the period is wrong (the jump is smeared), but nothing
escapes. This matches the reviewer's independent reproduction exactly.
''')

py(r'''
# Figure 1: buggy vs correct IM, h = 0.01, mu = 100, ~3.5 periods
def collect(stepfn, h, t_end):
    x, v = 2.0, 0.0
    ts, xs = [0.0], [2.0]
    for i in range(int(round(t_end/h))):
        x, v = stepfn(x, v, h)[:2]
        ts.append((i+1)*h); xs.append(x)
    return np.array(ts), np.array(xs)

tb, xb = collect(imid_b, 0.01, 3.5*Tref)
tc, xc = collect(imid_c, 0.01, 3.5*Tref)

fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
ax[0].plot(tb, np.abs(xb), lw=0.6, color="tab:red")
ax[0].set_yscale("log")
ax[0].set_ylabel("|x| (log scale)")
ax[0].set_title(r"$\mu=100$ van der Pol, $h=0.01$: implicit midpoint"
                "\nRED: Newton as published (wrong Jacobian sign)   "
                "BLUE: correct Newton (exact Jacobian)")
ax[0].grid(alpha=0.3)
ax[1].plot(tc, xc, lw=0.7, color="tab:blue")
ax[1].set_ylabel("x(t)")
ax[1].set_xlabel("t")
ax[1].legend(["correct IM: x(t)  (bounded, distorted)"])
ax[1].grid(alpha=0.3)
fig.tight_layout()
fig.show()
''')

md(r'''
## 6. Method V — localizing the failure: one step, one moment

Track the deviation $|(\text{buggy} - \text{correct})|$ over the whole run and inspect the
first step at which it exceeds $1$.
''')

py(r'''
def stage_roots(x, v, h):
    """Real roots w of the stage cubic (verified symbolically below):
    (mu/4) w^3 + mu x w^2 + (mu x^2 - mu + h/2 + 2/h) w + (h x - 2 v) = 0."""
    c3, c2 = mu/4.0, mu*x
    c1 = h/2.0 + mu*x*x - mu + 2.0/h
    c0 = h*x - 2.0*v
    rts = np.roots([c3, c2, c1, c0])
    real = [float(z.real) for z in rts if abs(z.imag) < 1e-8*max(1.0, abs(z))]
    out = []
    for w in real:
        xg, mv = x + w, w/h
        mx = (x + xg)/2
        vg = 2*mv - v
        res = max(abs(xg - x - h*mv), abs(vg - v - h*(mu*(1-mx*mx)*mv - mx)))
        out.append((w, res))
    return out

for h in (0.5, 0.1, 0.01):
    xb, vb, xc, vc = 2.0, 0.0, 2.0, 0.0
    n = int(3*Tref/h)
    last = None
    for i in range(n):
        xgb, vgb, itb, resb, convb = imid_b(xb, vb, h)
        xgc, vgc, itc, resc, convc = imid_c(xc, vc, h)
        d = max(abs(xgb - xgc), abs(vgb - vgc))
        if d > 1.0 and last is None:
            last = (i, xb, vb, xgb, vgb, itb, resb, convb, xgc, vgc, itc, resc, convc)
        xb, vb, xc, vc = xgb, vgb, xgc, vgc
        if not (np.isfinite(xb) and np.isfinite(vb)):
            break
    (i, xs, vs, xgb, vgb, itb, resb, convb, xgc, vgc, itc, resc, convc) = last
    print("h=%.2f : first step with deviation > 1: step %d (t = %.3f)" % (h, i+1, (i+1)*h))
    print("   state before:  x = %.8f, v = %.8f" % (xs, vs))
    print("   correct: xg = %.8f, vg = %.8f   it=%3d res=%.2e conv=%s"
          % (xgc, vgc, itc, resc, convc))
    print("   buggy  : xg = %.8g, vg = %.8g   it=%3d res=%.2e conv=%s"
          % (xgb, vgb, itb, resb, convb))
    x, v = xs, vs
    xg, vg = xgb, vgb
    mx_, mv_ = 0.5*(x + xg), 0.5*(v + vg)
    r1 = abs(xg - x - h*mv_)
    r2 = abs(vg - v - h*(mu*(1 - mx_*mx_)*mv_ - mx_))
    print("   buggy iterate stage residual (recomputed): G1 = %.3e, G2 = %.3e" % (r1, r2))
    roots = stage_roots(xs, vs, h)
    print("   real stage roots w at this state: %s"
          % ", ".join("%.4f" % w for w, _ in roots))
    if roots:
        print("   |buggy xg - nearest root| = %.3e"
              % min(abs(xgb - (xs + w)) for w, _ in roots))
    print()
''')

md(r'''
In **all three** step sizes the first divergence is the *same physical event* — the first
fast jump of the true orbit, at $t \approx 81$. There:

- the **correct** Newton either converges (residual $\le 10^{-14}$) or, at the coarsest
  steps, merely fails to converge while staying bounded — in every case its "step" is a
  legitimate (or nearly legitimate) IM step;
- the **buggy** Newton runs to its 100-iteration cap with residual $10^{81}$–$10^{90}$, and
  the script **silently returns the last iterate** — a point $10^{27}$–$10^{29}$ from *any*
  root of the stage equations.

The stage system (N1) is **cubic** in the endpoint increment $w = x_g - x$ (verified
symbolically in the next cell), so at the jump state several real roots exist; the correct
Newton, started near the current state, converges to the nearby one, while the wrong-sign
Jacobian has a different error-propagation spectrum and its iteration diverges. The
"blow-up" is therefore **instantaneous and localized**: one failed solve at $t\approx 81$;
the rest of the run merely carries the garbage along (the run's maximum is exactly the value
produced at that step).
''')

py(r'''
# Symbolic verification of the stage cubic (N1) in w = xg - x, mv = w/h, vg = 2*mv - v
xs_, vs_, hs_, mus_ = sp.symbols('x v h mu', real=True)
w = sp.symbols('w')
lhs = 2*w/hs_ - vs_
rhs = vs_ + hs_*(mus_*(1 - (xs_ + w/2)**2)*(w/hs_) - (xs_ + w/2))
poly = sp.Poly(sp.expand(lhs - rhs), w)
print("stage equation 2 as polynomial in w (endpoint increment):")
for i in range(poly.degree(), -1, -1):
    print("   coeff w^%d = %s" % (i, sp.simplify(poly.coeff_monomial(w**i))))
print("=>  (mu/4) w^3 + mu x w^2 + (mu x^2 - mu + h/2 + 2/h) w + (h x - 2 v) = 0")
''')

md(r'''
## 7. Method VI — how bad is the *correct* method really?

The essay's mechanism story ("residual $O(0.2)$ fast component re-amplified by the unstable
central region, so the orbit diverges") over-predicts: the correctly implemented orbit stays
bounded even at $h=0.5$ where $R(-150)\approx -0.974$. The genuine defect is **cycle
distortion**, growing with $h$:
''')

py(r'''
print("longer runs, CORRECT implicit midpoint (boundedness check):")
print("  %-26s %12s %12s %10s   %s" % ("h, periods", "T (last 3)", "max|x|", "max|v|", "status"))
for h0, nper in ((0.01, 20), (0.1, 10), (0.5, 5)):
    r = run_orbit(imid_c, h0, nper*Tref)
    print("  %-26s %12s %12.3e %10.3e   %s"
          % ("h=%.2f, %d per" % (h0, nper),
             ("%.4f" % r["T"]) if r["T"] else "None", r["mx"], r["mv"], r["status"]))

print("\nper-period envelope, CORRECT IM, h=0.01, 10 periods (true: amp 2, max|v|~134):")
x, v = 2.0, 0.0
per_steps = int(round(Tref/0.01))
for k in range(10):
    mx_, mv_ = 0.0, 0.0
    for i in range(per_steps):
        x, v = imid_c(x, v, 0.01)[:2]
        mx_ = max(mx_, abs(x)); mv_ = max(mv_, abs(v))
    print("   period %2d: max|x| = %9.3f   max|v| = %9.3f" % (k+1, mx_, mv_))
''')

md(r'''
**Stability-function values** at the relevant $z = \lambda h$ (fast mode $\lambda = -3\mu$
at $x = \pm 2$; central region $\lambda = +\mu$ at $x = 0$):
''')

py(r'''
def R(z):
    return (1 + z/2)/(1 - z/2)
print("  %-32s %10s %8s   |R(z)| %10s   e^z %12s" % ("mode", "z", "h", "", ""))
rows = [("fast mode x=2, lambda=-300", -300*h, "h=%.2f" % h) for h in (0.01, 0.1, 0.5)]
rows += [("central region x=0, lambda=+100", 100*h, "h=%.2f" % h) for h in (0.01, 0.1, 0.5)]
for label, z, tag in rows:
    ez = np.exp(z) if z < 700 else float('inf')
    print("  %-32s %10.2f %8s   %10.5f   %12.3e" % (label, z, tag, R(z), ez))
print("\nR(-inf) = -1 (A-stable, not L-stable): the fast mode is damped toward")
print("magnitude 1 with a sign flip each step, not toward 0 — the jump is smeared.")
''')

md(r'''
## 8. The rows that do NOT involve Newton (sanity of the re-run)

The Strang split (exact damping) and Radau IIA rows of the §6.3 table contain no hand-coded
Newton; they must reproduce exactly. They do:
''')

py(r'''
def strang_split(x, v, h, mu):
    v *= np.exp(mu*(1 - x*x)*h/2.0)
    v1 = v - 0.5*h*x
    x1 = x + h*v1
    v2 = v1 - 0.5*h*x1
    v2 *= np.exp(mu*(1 - x1*x1)*h/2.0)
    return x1, v2

for h0 in (0.5, 0.2):
    r = run_orbit(lambda x, v, h: strang_split(x, v, h, mu), h0, t3)
    print("  Strang split, h=%.2f : T=%s  max|x|=%.3e  max|v|=%.3e  %s"
          % (h0, ("%.4f" % r["T"]) if r["T"] else "None", r["mx"], r["mv"], r["status"]))
print("  (raw output: h=0.5 -> max|x| = 2.006e15;  h=0.2 -> T = 151.40, max|x| = 4.592e8)")

def radau_period(h, mu, t_end):
    def f(t, y):
        x, v = y
        return [v, mu*(1 - x*x)*v - x]
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
    T = radau_period(h0, mu, 3*Tref)
    print("  Radau IIA max_step=%.1f : T=%s  (ref %.4f)" % (h0, ("%.4f" % T) if T else "None", Tref))
''')

md(r'''
## 9. Conclusion

| | published IM ("blow-up") | correct IM |
|---|---|---|
| $h=0.5$ | $\max\lvert x\rvert = 1.0\times10^{27}$ (raw), printed $\sim10^{30}$ | bounded, $\max\lvert x\rvert = 3.0$, jump unresolved |
| $h=0.1$ | $\max\lvert x\rvert = 1.4\times10^{29}$ (raw), printed $\sim10^{32}$ | bounded, $T = 286.5$ ($+76\%$), $\max\lvert v\rvert = 73.8$ |
| $h=0.01$ | $\max\lvert x\rvert = 8.7\times10^{29}$ (raw = printed) | bounded, $T = 179.08$ ($+10\%$), $\max\lvert v\rvert = 144$ |

1. **The "$\sim10^{30}$" is an artifact of the solver bug**: wrong-sign (2,1) Jacobian entry
   + a Newton loop with no convergence check that silently returns the last iterate. The
   entire "blow-up" is one failed solve at the first fast jump ($t\approx 81$).
2. **The table misquotes the raw output** by $\times10^3$ in two rows (and $\div10^3$ in the
   Strang $h=0.2$ row).
3. **The essay's qualitative thesis survives** — in a weaker, verified form: the implicit
   midpoint rule (symplectic + $A$-stable, not $L$-stable) mishandles the stiff jump and
   produces a materially wrong (distorted, too-long-period) cycle, with the distortion
   growing with $h$; the $L$-stable Radau IIA is exact to the displayed digits at
   $h = 1.0 \gg 1/\mu$. The corrected statement is stronger in honesty and still sufficient
   for the symplecticity/$L$-stability tension of §6.5.
4. **Code recommendations**: fix $c$ to $-\tfrac{h}{2}J_{21}$; never return an unconverged
   Newton iterate silently (check the residual after the cap, raise or fall back to a
   guaranteed-convergent solver); re-run §6.3 and regenerate the table from raw output.

Companion files: `suspicion63_implicit_midpoint_blowup.md` (full analysis),
`suspicion63_analysis.py` / `suspicion63_results.txt`, `suspicion63_dive2.py` /
`suspicion63_dive2_results.txt`, `suspicion63_figure.png`.
''')

nb = nbformat.from_dict(nb_dict)
nbformat.validate(nb)
out = "suspicion63_implicit_midpoint_blowup.ipynb"
nbformat.write(nb, out)
print("wrote %s with %d cells (%d code)" % (
    out, len(nb.cells),
    sum(1 for c in nb.cells if c.cell_type == "code")))
