"""
Computation 3 (FIXED): implicit midpoint on stiff van der Pol.
Bug fix: Newton 2x2 update sign.  J=[[a,b],[c,d]], G=0:
    x_new = x - ( d*G1 - b*G2)/det
    y_new = y - (-c*G1 + a*G2)/det
Run: ./.venv/Scripts/python.exe geometric_integration_computation3.py
"""
import numpy as np
from scipy.integrate import solve_ivp

OUT=[]
def p(s=""): OUT.append(s); print(s)

def vdp_f(x,v,mu): return (v, mu*(1-x*x)*v - x)

def imid_step(x,v,h,f,mu):
    """Implicit midpoint (Gauss-Legendre 1-stage): symplectic & A-stable.
    Solve (x*,v*)=(x,v)+h f((x+x*)/2,(v+v*)/2) by Newton on G(x*,v*)=0."""
    xg,vg=x,v
    for _ in range(100):
        mx,mv=0.5*(x+xg),0.5*(v+vg)
        fx,fv=f(mx,mv)
        G1=xg-x-h*fx; G2=vg-v-h*fv
        if abs(G1)<1e-13 and abs(G2)<1e-13: return xg,vg
        J11,J12=0.0,1.0
        J21,J22=(-2.0*mu*mx*mv-1.0),(mu*(1-mx*mx))
        a=1.0-0.5*h*J11; b=-0.5*h*J12; c=0.5*h*J21; d=1.0-0.5*h*J22
        det=a*d-b*c
        xg-= ( d*G1 - b*G2)/det
        vg-= (-c*G1 + a*G2)/det
    return xg,vg

# sanity: harmonic oscillator implicit midpoint should stay on E=const (bounded)
def harm_f(x,v): return (v,-x)
x,v=1.0,0.0; E0=0.5; mxe=0.0
for i in range(100000):
    x,v=imid_step(x,v,0.1,harm_f,1.0)
    mxe=max(mxe,abs(0.5*x*x+0.5*v*v-E0))
p("SANITY harmonic IM h=0.1, 1e5 steps: max|E-E0| = %.3e (should be ~1e-5, bounded)" % mxe)

def ref_period(mu):
    def f(t,y):
        x,v=y; return [v, mu*(1-x*x)*v-x]
    def cross(t,y): return y[0]
    cross.terminal=False; cross.direction=1
    sol=solve_ivp(f,(0.0,4.0*mu*1.614+5.0),[2.0,0.0],method="Radau",
                  rtol=1e-10,atol=1e-10,events=cross,max_step=0.2)
    per=np.diff(sol.t_events[0])
    return float(np.mean(per[-3:]))

mu=100.0
Tref=ref_period(mu)
p("")
p("="*70)
p("Implicit midpoint on van der Pol, mu=100.  True: T=%.3f, amp->2, max|v|~134"%Tref)
p("Run IM for 400 time units; report (max|x|, max|v|) and final state (should be ~bounded).")
p("="*70)
p(" h      max|x|    max|v|   final x   final v")
for h in [0.5,0.1]:
    x,v=2.0,0.0; mx=abs(x); mv=abs(v)
    for i in range(int(400/h)):
        x,v=imid_step(x,v,h,lambda a,b: vdp_f(a,b,mu),mu)
        mx=max(mx,abs(x)); mv=max(mv,abs(v))
    p("%.3f  %9.4f  %9.4f  %9.4f  %9.4f"%(h,mx,mv,x,v))

def period_im(h,mu,t_end):
    x,v=2.0,0.0; xprev=x; cross=[]; mx=abs(x); mv=abs(v)
    for i in range(int(t_end/h)):
        x,v=imid_step(x,v,h,lambda a,b: vdp_f(a,b,mu),mu)
        mx=max(mx,abs(x)); mv=max(mv,abs(v))
        if xprev<0.0 and x>=0.0 and v>0.0: cross.append((i+1)*h)
        xprev=x
    if len(cross)<2: return None,mx,mv
    return float(np.mean(np.diff(cross)[-3:])),mx,mv

p("")
p("="*70)
p("Implicit midpoint period + envelope, mu=100 (stiff) and mu=1 (non-stiff).")
p("  mu=100 ref T=%.3f (amp 2, max|v|~134);  mu=1 ref T=6.6633 (amp~2.008)"%Tref)
p("="*70)
p(" mu=100:")
for h in [0.5,0.1]:
    T,mx,mv=period_im(h,mu,3.0*Tref)
    p("   h=%.2f  T=%s  max|x|=%.4f  max|v|=%.3f"%(h,("%.4f"%T) if T else "None",mx,mv))
p(" mu=1:")
for h in [0.01,0.005,0.001]:
    T,mx,mv=period_im(h,1.0,30.0)
    p("   h=%.4f  T=%s  max|x|=%.4f  max|v|=%.4f"%(h,("%.6f"%T) if T else "None",mx,mv))

with open("geometric_integration_results3.txt","w") as fh: fh.write("\n".join(OUT)+"\n")
p("\n[written]")
