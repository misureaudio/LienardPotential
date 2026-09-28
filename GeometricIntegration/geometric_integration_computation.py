"""
Geometric integration vs. the Lienard equation: grounding computations.

Every numerical value cited in geometric_integration_lienard_essay_v0.md is
produced here. Run with the project venv:
    ./.venv/Scripts/python.exe geometric_integration_computation.py

Computations:
  (D) Harmonic oscillator: Storer-Verlet (symplectic) vs RK4 energy error, long time.
  (A) Duffing double-well (conservative): symplectic vs RK4 energy drift.
  (B1) Van der Pol period, mu=1: RK4 vs symplectic-dissipative Strang, h-refinement.
  (B3) Van der Pol relaxation period (implicit Radau): T(mu) ~ mu(3 - 2 ln 2).
  (E) Van der Pol discrete Poincare map (mu=0.1): fixed point vs h (discrete limit cycle).
  (S) Symbolic: Storer-Verlet map is symplectic (det J = 1, and M^T J M = J).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import sympy as sp

OUT = []
def p(s=""):
    OUT.append(s)
    print(s)

# ----------------------------------------------------------------------
# Basic one-step methods
# ----------------------------------------------------------------------
def sv_step(x, v, h, Vp):
    """Storer-Verlet (leapfrog) for separable H = v^2/2 + V(x). Symplectic, 2nd order."""
    v1 = v - 0.5 * h * Vp(x)
    x1 = x + h * v1
    v2 = v1 - 0.5 * h * Vp(x1)
    return x1, v2

def rk4_step(x, v, h, f):
    """Classic RK4 for (x', v') = f(x, v). 4th order, not symplectic."""
    k1x, k1v = f(x, v)
    k2x, k2v = f(x + 0.5*h*k1x, v + 0.5*h*k1v)
    k3x, k3v = f(x + 0.5*h*k2x, v + 0.5*h*k2v)
    k4x, k4v = f(x + h*k3x, v + h*k3v)
    x1 = x + (h/6.0)*(k1x + 2*k2x + 2*k3x + k4x)
    v1 = v + (h/6.0)*(k1v + 2*k2v + 2*k3v + k4v)
    return x1, v1

# Systems
def harm_f(x, v):   return (v, -x)
def duff_f(x, v, a, b): return (v, -(a*x + b*x**3))
def vdp_f(x, v, mu):    return (v, mu*(1 - x*x)*v - x)

def strang_vdp(x, v, h, mu):
    """Symplectic-dissipative Strang splitting for van der Pol.
    Hamiltonian part H = x^2/2 + v^2/2 (Vp = x); damping dv/dt = mu(1-x^2) v (exact exp)."""
    v = v * np.exp(mu*(1 - x*x) * h/2.0)      # half damping (exact)
    v1 = v - 0.5*h*x
    x1 = x + h*v1
    v2 = v1 - 0.5*h*x1
    v2 = v2 * np.exp(mu*(1 - x1*x1) * h/2.0)  # half damping (exact)
    return x1, v2

# ----------------------------------------------------------------------
# (D) Harmonic oscillator energy error
# ----------------------------------------------------------------------
def energy_harm(x, v): return 0.5*x*x + 0.5*v*v

def run_harmonic(h, t_end, method):
    x, v = 1.0, 0.0
    E0 = energy_harm(x, v)
    n = int(round(t_end/h))
    max_err = 0.0
    for i in range(n):
        if method == "sv":
            x, v = sv_step(x, v, h, lambda xx: xx)
        else:
            x, v = rk4_step(x, v, h, harm_f)
        e = abs(energy_harm(x, v) - E0)
        if e > max_err: max_err = e
    return max_err, energy_harm(x, v) - E0

p("="*70)
p("(D) Harmonic oscillator  x'' + x = 0,  x(0)=1, v(0)=0,  E0 = 1/2")
p("    t_end = 1000 (~159 periods).  max|E-E0| over the run.")
p("="*70)
p(f"{'h':>8} {'method':>14} {'max|E-E0|':>14} {'E-E0 at tend':>14}")
for h in [0.1, 0.01]:
    for m, name in [("sv","Storer-Verlet"), ("rk4","RK4")]:
        mx, end = run_harmonic(h, 1000.0, m)
        p(f"{h:>8.3f} {name:>14} {mx:>14.6e} {end:>+14.6e}")

# ----------------------------------------------------------------------
# (A) Duffing double-well energy drift (conservative)
# ----------------------------------------------------------------------
a_d, b_d = -1.0, 1.0   # V = -x^2/2 + x^4/4, wells at x = +-1
def energy_duff(x, v): return 0.5*a_d*x*x + 0.5*v*v + 0.25*b_d*x**4

def run_duffing(h, t_end, method):
    x, v = 0.5, 0.0    # in the right well (E0 < 0 = saddle level)
    E0 = energy_duff(x, v)
    n = int(round(t_end/h))
    max_err = 0.0
    ts, es = [], []
    for i in range(n):
        if method == "sv":
            x, v = sv_step(x, v, h, lambda xx: a_d*xx + b_d*xx**3)
        else:
            x, v = rk4_step(x, v, h, lambda xx, vv: duff_f(xx, vv, a_d, b_d))
        t = (i+1)*h
        e = energy_duff(x, v) - E0
        ae = abs(e)
        if ae > max_err: max_err = ae
        if i % max(1, n//200) == 0:
            ts.append(t); es.append(e)
    ts = np.array(ts); es = np.array(es)
    slope = np.polyfit(ts, es, 1)[0]   # secular drift rate
    return max_err, es[-1], slope, E0

p("")
p("="*70)
p("(A) Duffing double-well  x'' - x + x^3 = 0  (a=-1, b=1),  x(0)=0.5, v(0)=0")
p("    V = -x^2/2 + x^4/4 (wells at x=+-1).  t_end = 2000 (~840 periods).")
p("    E0 = %.8f (below saddle level 0 => closed orbit in right well)." % run_duffing(0.05, 1.0, "sv")[3])
p("    max|E-E0|, E-E0 at tend, and secular drift rate (linear fit of E-E0 vs t).")
p("="*70)
p(f"{'h':>8} {'method':>14} {'max|E-E0|':>14} {'E-E0@tend':>14} {'drift rate':>14}")
for h in [0.05, 0.01]:
    for m, name in [("sv","Storer-Verlet"), ("rk4","RK4")]:
        mx, end, slope, _ = run_duffing(h, 2000.0, m)
        p(f"{h:>8.3f} {name:>14} {mx:>14.6e} {end:>+14.6e} {slope:>+14.6e}")

# ----------------------------------------------------------------------
# (B1) Van der Pol period, mu = 1
# ----------------------------------------------------------------------
def period_vdp(h, mu, method, t_span=(0.0, 60.0)):
    """Integrate van der Pol and measure the period from successive
    upward crossings of x = 0 (v > 0). Returns list of periods."""
    x, v = 2.0, 0.0
    t = 0.0
    crossings = []
    x_prev, v_prev = x, v
    n = int(round((t_span[1]-t_span[0])/h))
    for i in range(n):
        if method == "rk4":
            x, v = rk4_step(x, v, h, lambda xx, vv: vdp_f(xx, vv, mu))
        else:
            x, v = strang_vdp(x, v, h, mu)
        t += h
        # upward crossing of x=0 with v>0
        if x_prev < 0.0 and x >= 0.0 and v > 0.0:
            frac = -x_prev/(x - x_prev)
            t_c = (t-h) + frac*h
            crossings.append(t_c)
        x_prev, v_prev = x, v
    periods = np.diff(crossings)
    return periods

def reference_period(mu):
    """High-accuracy period via implicit Radau with event detection."""
    def f(t, y):
        x, v = y
        return [v, mu*(1 - x*x)*v - x]
    def cross(t, y): return y[0]
    cross.terminal = False
    cross.direction = 1
    sol = solve_ivp(f, (0.0, 60.0), [2.0, 0.0], method="Radau",
                    rtol=1e-11, atol=1e-11, events=cross, max_step=0.5)
    per = np.diff(sol.t_events[0])
    return float(np.mean(per[1:]))  # discard first (possible transient)

p("")
p("="*70)
p("(B1) Van der Pol period, mu = 1.  Reference (Radau, rtol=1e-11):")
ref1 = reference_period(1.0)
p("    T_ref = %.8f" % ref1)
p("    Period from successive upward x=0 crossings (mean of last 5).")
p("="*70)
p(f"{'h':>9} {'RK4 period':>14} {'Strang period':>16} {'|RK4-ref|':>12} {'|Str-ref|':>12}")
for h in [0.01, 0.005, 0.001]:
    pr = period_vdp(h, 1.0, "rk4")
    ps = period_vdp(h, 1.0, "strang")
    mr = float(np.mean(pr[-5:])); ms = float(np.mean(ps[-5:]))
    p(f"{h:>9.4f} {mr:>14.8f} {ms:>16.8f} {abs(mr-ref1):>12.3e} {abs(ms-ref1):>12.3e}")

# ----------------------------------------------------------------------
# (B3) Van der Pol relaxation period (implicit)
# ----------------------------------------------------------------------
def relaxation_period(mu):
    def f(t, y):
        x, v = y
        return [v, mu*(1 - x*x)*v - x]
    def cross(t, y): return y[0]
    cross.terminal = False
    cross.direction = 1
    t_end = 4.0 * mu * 1.614 + 5.0   # ~4 periods
    sol = solve_ivp(f, (0.0, t_end), [2.0, 0.0], method="Radau",
                    rtol=1e-10, atol=1e-10, events=cross, max_step=0.2)
    per = np.diff(sol.t_events[0])
    return float(np.mean(per[-3:]))

C = 3.0 - 2.0*np.log(2.0)
p("")
p("="*70)
p("(B3) Van der Pol relaxation period (Radau, rtol=1e-10), mu large.")
p("    Asymptote: T(mu) = mu*(3 - 2 ln 2) + O(1),   3 - 2 ln 2 = %.8f" % C)
p("="*70)
p(f"{'mu':>6} {'T (measured)':>14} {'T/mu':>10} {'T/(mu*C)':>12}")
for mu in [10.0, 50.0, 100.0]:
    T = relaxation_period(mu)
    p(f"{mu:>6.0f} {T:>14.4f} {T/mu:>10.5f} {T/(mu*C):>12.5f}")

# ----------------------------------------------------------------------
# (E) Van der Pol discrete Poincare map (mu = 0.1)
# ----------------------------------------------------------------------
def pmap_v(v0, mu, h, method):
    """Poincare map: start at (x=0, v=v0), return v at next upward x=0 crossing."""
    x, v = 0.0, v0
    x_prev, v_prev = x, v
    n = int(round(12.0/h))   # one period ~ 2*pi ~ 6.28; allow headroom
    for i in range(n):
        if method == "rk4":
            x, v = rk4_step(x, v, h, lambda xx, vv: vdp_f(xx, vv, mu))
        else:
            x, v = strang_vdp(x, v, h, mu)
        if x_prev < 0.0 and x >= 0.0 and v > 0.0:
            frac = -x_prev/(x - x_prev)
            return v_prev + frac*(v - v_prev)
        x_prev, v_prev = x, v
    return None

p("")
p("="*70)
p("(E) Van der Pol discrete Poincare map, mu = 0.1 (Strang splitting).")
p("    Fixed point v* of P(v)=v = discrete limit cycle (x=0 section, v>0).")
p("    Continuous limit: v* -> 2 (circle of radius 2), amplitude ~2.0001.")
p("    g(v) = P(v) - v; v* found by brentq on [1.5, 2.5].")
p("="*70)
mu_e = 0.1
p(f"{'h':>9} {'v* (fixed pt)':>16} {'|v* - 2|':>12}")
for h in [0.01, 0.005, 0.001, 0.0005]:
    g = lambda v: pmap_v(v, mu_e, h, "strang") - v
    try:
        vs = brentq(g, 1.5, 2.5, xtol=1e-12)
        p(f"{h:>9.4f} {vs:>16.10f} {abs(vs-2.0):>12.3e}")
    except Exception as ex:
        p(f"{h:>9.4f}   (brentq failed: {ex})")

# ----------------------------------------------------------------------
# (S) Symbolic: Storer-Verlet is symplectic
# ----------------------------------------------------------------------
p("")
p("="*70)
p("(S) Symbolic check: Storer-Verlet map for the harmonic oscillator is symplectic.")
p("    In 2D, a map M is symplectic (M^T J M = J, J=[[0,1],[-1,0]]) iff det M = 1.")
p("="*70)
x, v, h = sp.symbols('x v h')
v1 = v - sp.Rational(1,2)*h*x
x1 = x + h*v1
v2 = v1 - sp.Rational(1,2)*h*x1
x1 = sp.expand(x1); v2 = sp.expand(v2)
M = sp.Matrix([[sp.diff(x1, x), sp.diff(x1, v)],
               [sp.diff(v2, x), sp.diff(v2, v)]])
J = sp.Matrix([[0, 1], [-1, 0]])
detM = sp.simplify(M.det())
sympl = sp.simplify(M.T*J*M - J)
p("    x1 = %s" % sp.sstr(x1))
p("    v2 = %s" % sp.sstr(v2))
p("    det(M) = %s" % detM)
p("    M^T J M - J = %s   (zero => symplectic)" % sp.sstr(sympl))

# ----------------------------------------------------------------------
with open("geometric_integration_results.txt", "w") as fh:
    fh.write("\n".join(OUT) + "\n")
p("")
p("[results written to geometric_integration_results.txt]")
