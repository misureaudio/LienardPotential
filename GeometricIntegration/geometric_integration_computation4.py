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

# --- Implicit midpoint (symplectic, A-stable, not L-stable) ---
def imid_step(x,v,h,mu):
    xg,vg=x,v
    for _ in range(100):
        mx,mv=0.5*(x+xg),0.5*(v+vg)
        fx,fv=vdp_f(mx,mv,mu)
        G1=xg-x-h*fx; G2=vg-v-h*fv
        if abs(G1)<1e-13 and abs(G2)<1e-13: return xg,vg
        J11,J12=0.0,1.0
        J21,J22=(-2.0*mu*mx*mv-1.0),(mu*(1-mx*mx))
        a=1.0-0.5*h*J11; b=-0.5*h*J12; c=0.5*h*J21; d=1.0-0.5*h*J22
        det=a*d-b*c
        xg-= ( d*G1 - b*G2)/det
        vg-= (-c*G1 + a*G2)/det
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
    x,v=2.0,0.0; xprev=x; cross=[]; mx=abs(x); mv=abs(v)
    for i in range(int(t_end/h)):
        x,v=stepfn(x,v,h,mu)
        if not (np.isfinite(x) and np.isfinite(v)): return None,mx,mv,i+1,"DIVERGED"
        mx=max(mx,abs(x)); mv=max(mv,abs(v))
        if xprev<0.0 and x>=0.0 and v>0.0: cross.append((i+1)*h)
        xprev=x
    if len(cross)<2: return None,mx,mv,int(t_end/h),"no-crossings"
    return float(np.mean(np.diff(cross)[-3:])),mx,mv,int(t_end/h),"ok"

mu=100.0
Tref=ref_period(mu)
p("="*72)
p("Stiff van der Pol, mu=100.  True: T=%.3f, amplitude->2, max|v|~134"%Tref)
p("All methods start at (x,v)=(2,0), run 3 periods (~%d time units)."% (3*Tref))
p("="*72)
p("%-38s %10s %10s %12s %12s"%("method (h)","T","max|x|","max|v|","status"))
rows=[]
rows.append(("Implicit midpoint, h=0.5", run(imid_step,0.5,mu,3*Tref)))
rows.append(("Implicit midpoint, h=0.1", run(imid_step,0.1,mu,3*Tref)))
rows.append(("Implicit midpoint, h=0.01", run(imid_step,0.01,mu,3*Tref)))
rows.append(("Strang split (exact damp), h=0.5", run(strang_split,0.5,mu,3*Tref)))
rows.append(("Strang split (exact damp), h=0.2", run(strang_split,0.2,mu,3*Tref)))
for name,(T,mx,mv,n,st) in rows:
    p("%-38s %10s %10.4f %12.3f %12s"%(name, ("%.4f"%T) if T else "None", mx, mv, st))

# Radau IIA via scipy (L-stable, not symplectic) as the benchmark that WORKS
p("")
p("Benchmark that works in the stiff regime (L-stable, NOT symplectic):")
def radau_period(h,mu,t_end):
    def f(t,y):
        x,v=y; return [v, mu*(1-x*x)*v-x]
    def cross(t,y): return y[0]
    cross.terminal=False; cross.direction=1
    sol=solve_ivp(f,(0.0,t_end),[2.0,0.0],method="Radau",rtol=1e-9,atol=1e-9,
                  events=cross,max_step=h)
    if len(sol.t_events[0])<2: return None
    return float(np.mean(np.diff(sol.t_events[0])[-3:]))
for h in [0.5,1.0]:
    T=radau_period(h,mu,3*Tref)
    p("  Radau IIA, max_step=%s :  T=%s   (ref %.3f)"%(h,("%.4f"%T) if T else "None",Tref))

p("")
p("Reference: Tref=%.4f, |T_Radau - Tref| ~ %s"%(Tref, "%.2e"%(abs(radau_period(1.0,mu,3*Tref)-Tref)) if radau_period(1.0,mu,3*Tref) else "n/a"))
with open("geometric_integration_results4.txt","w") as fh: fh.write("\n".join(OUT)+"\n")
p("\n[written]")
