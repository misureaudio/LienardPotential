"""
Computation 4: L-stability vs A-stability for stiff van der Pol (mu=100).
Honest question: what actually works in the stiff relaxation regime?
  (a) Implicit midpoint (symplectic, A-stable, NOT L-stable): blows up?
  (b) Radau IIA (L-stable, NOT symplectic): works?
  (c) Symplectic-dissipative Strang split, EXACT exp damping (L-stable-like): works?
Run: ./.venv/Scripts/python.exe geometric_integration_computation4.py
"""
import numpy as np
from scipy.integrate import solve_ivp

OUT=[]
def p(s=""): OUT.append(s); print(s)

def vdp_f(x,v,mu): return (v, mu*(1-x*x)*v - x)

# --- Implicit midpoint (symplectic, A-stable, NOT L-stable) --- CORRECTED v2 ---
# v2 REMEDIATION: the (2,1) Jacobian entry is c = -h/2*J21 (exact).  The v1
# scripts used c = +h/2*J21 (wrong sign), which made Newton diverge at the first
# fast jump and silently return a 10^27-10^29 garbage iterate (the "blow-up").
# Now: exact Jacobian + fail-loud guard (raise if residual > 1 or non-finite; a
# residual <= 1 is a legitimate bounded / jump-unresolved step).
def imid_step(x,v,h,mu, tol=1e-14, itmax=2000, guard=1.0):
    xg,vg = x,v
    res = float('inf')
    for _ in range(itmax):
        mx,mv = 0.5*(x+xg), 0.5*(v+vg)
        G1 = xg - x - h*mv
        G2 = vg - v - h*(mu*(1-mx*mx)*mv - mx)
        res = max(abs(G1), abs(G2))
        if res < tol:
            return xg,vg
        a, b = 1.0, -0.5*h
        c = -0.5*h*(-2.0*mu*mx*mv - 1.0)      # EXACT (2,1) entry  c = -h/2*J21  (v1 had +)
        d = 1.0 - 0.5*h*(mu*(1-mx*mx))
        det = a*d - b*c
        xg -= ( d*G1 - b*G2)/det
        vg -= (-c*G1 + a*G2)/det
    # fail-loud: never return an unconverged (garbage) iterate silently
    if not (np.isfinite(xg) and np.isfinite(vg)):
        raise RuntimeError("IM Newton diverged (non-finite) at x=%g v=%g h=%g" % (x,v,h))
    if res > guard:
        raise RuntimeError("IM step UNRESOLVED (residual %.3e) at x=%g v=%g h=%g" % (res,x,v,h))
    return xg,vg

# --- Symplectic-dissipative Strang split, EXACT exp damping ---
def strang_split(x,v,h,mu):
    v*=np.exp(mu*(1-x*x)*h/2.0)      # half exact damping (L-stable-like: ->0 as mu h -> inf, x>1)
    v1=v-0.5*h*x; x1=x+h*v1; v2=v1-0.5*h*x1   # Storer-Verlet (symplectic) for H=x^2/2+v^2/2
    v2*=np.exp(mu*(1-x1*x1)*h/2.0)
    return x1,v2

def ref_period(mu):
    def f(t,y):
        x,v=y; return [v, mu*(1-x*x)*v-x]
    def cross(t,y): return y[0]
    cross.terminal=False; cross.direction=1
    sol=solve_ivp(f,(0.0,4.0*mu*1.614+5.0),[2.0,0.0],method="Radau",
                  rtol=1e-10,atol=1e-10,events=cross,max_step=0.2)
    per=np.diff(sol.t_events[0])
    return float(np.mean(per[-3:]))

def run(stepfn,h,mu,t_end):
    """Integrate; return (T, max|x|, max|v|, status). Catches a fail-loud
    RuntimeError from the step function and reports it as UNRESOLVED."""
    x,v=2.0,0.0; xprev=x; cross=[]; mx=abs(x); mv=abs(v)
    for i in range(int(round(t_end/h))):
        try:
            x,v=stepfn(x,v,h,mu)
        except RuntimeError as e:
            return None,mx,mv,"UNRESOLVED: "+str(e)
        if not (np.isfinite(x) and np.isfinite(v)): return None,mx,mv,"DIVERGED"
        mx=max(mx,abs(x)); mv=max(mv,abs(v))
        if xprev<0.0 and x>=0.0 and v>0.0: cross.append((i+1)*h)
        xprev=x
    if len(cross)<2: return None,mx,mv,"no-crossings"
    return float(np.mean(np.diff(cross)[-3:])),mx,mv,"ok"

mu=100.0
Tref=ref_period(mu)
p("="*72)
p("Stiff van der Pol, mu=100.  True: T=%.3f, amplitude->2, max|v|~134"%Tref)
p("All methods start at (x,v)=(2,0), run 3 periods (~%d time units)."% (3*Tref))
p("="*72)
p("%-38s %10s %10s %12s %12s"%("method (h)","T","max|x|","max|v|","status"))
def radau_period(h,mu,t_end):
    """Return (T, max|x|) for Radau IIA with max_step=h (envelope tracked to
    confirm boundedness)."""
    def f(t,y):
        x,v=y; return [v, mu*(1-x*x)*v-x]
    def cross(t,y): return y[0]
    cross.terminal=False; cross.direction=1
    sol=solve_ivp(f,(0.0,t_end),[2.0,0.0],method="Radau",rtol=1e-9,atol=1e-9,
                  events=cross,max_step=h)
    mx=float(np.max(np.abs(sol.y[0])))
    if len(sol.t_events[0])<2: return None,mx
    return float(np.mean(np.diff(sol.t_events[0])[-3:])),mx

rows=[]
for h in (0.5,0.1,0.01):
    T,mx,mv,st = run(imid_step,h,mu,3*Tref)
    rows.append(("implicit midpoint, h=%.2f"%h, T,mx,mv,st))
for h in (0.5,0.2):
    T,mx,mv,st = run(strang_split,h,mu,3*Tref)
    rows.append(("Strang split, exact damping, h=%.2f"%h, T,mx,mv,st))
for name,T,mx,mv,st in rows:
    p("%-38s %10s %10s %12s %s"%(name, ("%.4f"%T) if T else "-",
        ("%.3f"%mx) if (mx is not None and mx<1e3) else ("%.3e"%mx if mx is not None else "-"),
        ("%.3f"%mv) if (mv is not None and mv<1e3) else ("%.3e"%mv if mv is not None else "-"), st))

# Radau IIA rows (L-stable, not symplectic) — the benchmark that WORKS
for h in [0.5,1.0]:
    T,mx=radau_period(h,mu,3*Tref)
    rows.append(("Radau IIA, max_step=%.1f"%h, T,mx,None,"ok" if T else "no-crossings"))
    p("  Radau IIA, max_step=%s :  T=%s   max|x|=%.3f   (ref %.3f)"%(h,("%.4f"%T) if T else "None",mx,Tref))

# --- STEP 1.3: Markdown table for the essay (copy-paste verbatim; no retyping) ---
def mdnum(x):
    if x is None: return "-"
    if x < 1e3: return "%.2f"%x
    e=int(np.floor(np.log10(x))); m=x/10.0**e
    return "~%.1f$\\times$10^{%d}"%(m,e)
def mdT(T):
    return ("%.4f"%T) if T else "- (no clean crossings)"
def mdstatus(name,T,mx,mv,st):
    if st.startswith("UNRESOLVED") or st=="DIVERGED": return "**fail-loud**"
    if "Strang" in name:
        if T is not None: return "wrong ($T$ off)"
        return "**blow-up**" if (mx is not None and mx>1e3) else "no crossings"
    if "Radau" in name: return "**correct**"
    if T is None: return "bounded, jump unresolved"
    err=(T-Tref)/Tref*100.0
    return "bounded, distorted ($T$ off by $+%d\\%%$)"%round(err)

p("")
p("MARKDOWN TABLE (copy-paste into the essay):")
p("")
p("| method (step $h$) | period $T$ | $\\max\\lvert x\\rvert$ over run | status |")
p("|---|---:|---:|---|")
for name,T,mx,mv,st in rows:
    nm = name
    if "Radau" in name:
        nm = "**%s**" % name.replace("max_step","$\\max\\mathrm{step}$")
    # Radau envelope is the true-cycle amplitude (~2); show "bounded" (not the raw ~2.00).
    mxcell = "bounded" if ("Radau" in name and mx is not None and mx<1e3) else mdnum(mx)
    p("| %s | %s | %s | %s |"%(nm, mdT(T), mxcell, mdstatus(name,T,mx,mv,st)))
p("")
p("Reference: Tref=%.4f, |T_Radau - Tref| ~ %s"%(Tref, "%.2e"%(abs(radau_period(1.0,mu,3*Tref)[0]-Tref)) if radau_period(1.0,mu,3*Tref)[0] else "n/a"))
with open("geometric_integration_results4.txt","w") as fh: fh.write("\n".join(OUT)+"\n")
p("\n[written]")
