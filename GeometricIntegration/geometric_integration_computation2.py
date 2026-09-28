"""
Computation 2: stiff regime + implicit midpoint (symplectic & A-stable).
Run: ./.venv/Scripts/python.exe geometric_integration_computation2.py
"""
import numpy as np
from scipy.integrate import solve_ivp
import sympy as sp

OUT = []
def p(s=""):
    OUT.append(s); print(s)

def vdp_f(x, v, mu): return (v, mu*(1 - x*x)*v - x)

def rk4_step(x, v, h, f):
    k1x, k1v = f(x, v)
    k2x, k2v = f(x + 0.5*h*k1x, v + 0.5*h*k1v)
    k3x, k3v = f(x + 0.5*h*k2x, v + 0.5*h*k2v)
    k4x, k4v = f(x + h*k3x, v + h*k3v)
    return (x + (h/6.0)*(k1x+2*k2x+2*k3x+k4x),
            v + (h/6.0)*(k1v+2*k2v+2*k3v+k4v))

def imid_step(x, v, h, f, mu):
    """Implicit midpoint (Gauss-Legendre 1-stage): symplectic & A-stable.
    Solve (x*,v*) = (x,v) + h f((x+x*)/2,(v+v*)/2) by Newton."""
    xg, vg = x, v
    for _ in range(50):
        mx, mv = 0.5*(x+xg), 0.5*(v+vg)
        fx, fv = f(mx, mv)
        rx = x + h*fx - xg
        rv = v + h*fv - vg
        if abs(rx) < 1e-13 and abs(rv) < 1e-13:
            return xg, vg
        # Jacobian of G(xg,vg)=(xg - x - h f(mid), vg - v - h f(mid))
        # d f(mid)/d(xg,vg) = 0.5 J_f(mid)
        jfx, jfv_x = fx, 0.0
        # f=(v, mu(1-x^2)v - x): J = [[0,1],[ -2 mu x v - 1, mu(1-x^2) ]]
        J11, J12 = 0.0, 1.0
        J21, J22 = (-2.0*mu*mx*mv - 1.0), (mu*(1 - mx*mx))
        a = 1.0 - 0.5*h*J11
        b = -0.5*h*J12
        c = -0.5*h*J21
        d = 1.0 - 0.5*h*J22
        det = a*d - b*c
        xg = xg - (d*rx - b*rv)/det
        vg = vg - (c*rx - a*rv)/det
    return xg, vg

def period_with(h, mu, stepfn, t_end, max_it=None):
    x, v = 2.0, 0.0
    x_prev = x
    cross = []
    n = int(round(t_end/h))
    for i in range(n):
        x, v = stepfn(x, v, h)
        if not (np.isfinite(x) and np.isfinite(v)):
            return None, i+1, "DIVERGED"
        if x_prev < 0.0 and x >= 0.0 and v > 0.0:
            cross.append((i+1)*h)
        x_prev = x
        if max_it and i >= max_it:
            break
    if len(cross) < 2:
        return None, n, "no crossings"
    per = np.diff(cross)
    return float(np.mean(per[-3:])), n, "ok"

# ----------------------------------------------------------------------
p("="*70)
p("(F1) Implicit midpoint: symplectic + A-stable.")
p("    One-step map for harmonic oscillator (x',v')=(v,-x):")
p("="*70)
x, v, h = sp.symbols('x v h')
# implicit midpoint: (x*,v*) = (x,v) + h*( (v+v*)/2 , -(x+x*)/2 )
xs, vs = sp.symbols('xs vs')
eq1 = sp.Eq(xs, x + h*(v+vs)/2)
eq2 = sp.Eq(vs, v + h*(-(x+xs)/2))
sol = sp.solve([eq1, eq2], [xs, vs], dict=True)[0]
xs_s = sp.simplify(sol[xs]); vs_s = sp.simplify(sol[vs])
M = sp.Matrix([[sp.diff(xs_s, x), sp.diff(xs_s, v)],
               [sp.diff(vs_s, x), sp.diff(vs_s, v)]])
J = sp.Matrix([[0,1],[-1,0]])
p("    x* = %s" % sp.sstr(xs_s))
p("    v* = %s" % sp.sstr(vs_s))
p("    det(M) = %s   (1 => symplectic)" % sp.simplify(M.det()))
p("    M^T J M - J = %s" % sp.sstr(sp.simplify(M.T*J*M - J)))
z = sp.symbols('z')
R = (1 + z/2)/(1 - z/2)
w = sp.symbols('w', real=True)
absR2 = sp.simplify(R.subs(z, sp.I*w).conjugate() * R.subs(z, sp.I*w))
p("    R(z)=(1+z/2)/(1-z/2);  |R(iw)|^2 = %s  (=1 => A-stable, |R|=1 on imag axis)" % absR2)
p("    R(-inf) = -1  =>  A-stable but NOT L-stable (Radau IIA is L-stable).")

# ----------------------------------------------------------------------
p("")
p("="*70)
p("(F2) Van der Pol, mu = 100 (stiff relaxation regime).  Reference T ~ 162.84.")
p("    RK4 stable step: h < 0.93/mu = 0.0093.  Compare RK4 h=0.005 vs IM h=0.5.")
p("="*70)
mu = 100.0
refT = 162.8371
t_end = 3.0*refT
# RK4 (expensive: ~ t_end/0.005 steps)
T_rk4, n_rk4, st = period_with(0.005, mu, lambda x,v,h: rk4_step(x,v,h,lambda a,b: vdp_f(a,b,mu)), t_end)
p("    RK4  h=0.005:  T = %s   (%s steps, %s)" %
  ("%.4f"%T_rk4 if T_rk4 else "None", n_rk4, st))
# Implicit midpoint (cheap: ~ t_end/0.5 steps)
T_im, n_im, st = period_with(0.5, mu, lambda x,v,h: imid_step(x,v,h,lambda a,b: vdp_f(a,b,mu), mu), t_end)
p("    IM   h=0.5  :  T = %s   (%s steps, %s)" %
  ("%.4f"%T_im if T_im else "None", n_im, st))
# RK4 with an UNSTABLE step h=0.5 (should blow up or be wrong)
T_bad, n_bad, st = period_with(0.5, mu, lambda x,v,h: rk4_step(x,v,h,lambda a,b: vdp_f(a,b,mu)), t_end, max_it=200000)
p("    RK4  h=0.5  :  T = %s   (%s steps, %s)   [h=0.5 > 0.0093 => unstable]" %
  ("%.4f"%T_bad if T_bad else "None", n_bad, st))
if T_rk4 and T_im:
    p("    steps/period:  RK4 ~ %.0f   vs   IM ~ %.0f   (ratio ~ %.0fx)" %
      (refT/0.005, refT/0.5, (refT/0.005)/(refT/0.5)))
    p("    period error:  RK4 |T-ref| = %.3e   IM |T-ref| = %.3e" %
      (abs(T_rk4-refT), abs(T_im-refT)))

# ----------------------------------------------------------------------
p("")
p("="*70)
p("(F3) Implicit midpoint orbit on mu=100 cycle: slow-manifold tracking.")
p("    Slow manifold y_s(x) = x/(mu(1-x^2)) for |x|>1.  Max |v - y_s(x)| on")
p("    the right slow branch (x in [1.05, 1.95]).")
p("="*70)
x, v = 2.0, 0.0
dev = []
n = int(round(2.0*refT/0.5))
for i in range(n):
    x, v = imid_step(x, v, 0.5, lambda a,b: vdp_f(a,b,mu), mu)
    if 1.05 < x < 1.95:
        ys = x/(mu*(1 - x*x))
        dev.append(abs(v - ys))
if dev:
    p("    max |v - y_s(x)| on right slow branch = %.3e   (IM, h=0.5)" % max(dev))

with open("geometric_integration_results2.txt","w") as fh:
    fh.write("\n".join(OUT)+"\n")
p("")
p("[written to geometric_integration_results2.txt]")
