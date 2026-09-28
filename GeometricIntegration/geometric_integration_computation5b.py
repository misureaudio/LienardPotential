"""
Computation 5b: resolve physical vs. Lienard-plane picture of the van der Pol
cycle at mu=0.1, and where the discrete Poincare fixed point lands.
Run: ./.venv/Scripts/python.exe geometric_integration_computation5b.py
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def F(x,mu): return mu*(x**3/3 - x)
def vdp(t,y,mu):
    x,v=y; return [v, mu*(1-x*x)*v - x]

def true_cycle(mu, rtol=1e-12):
    def cross(t,y,*a): return y[0]
    cross.terminal=False; cross.direction=1
    sol = solve_ivp(vdp,(0.0,20.0),[2.0,0.0],method="Radau",args=(mu,),
                    rtol=rtol,atol=rtol,events=cross,max_step=0.05)
    ye=sol.y_events[0]; vs=ye[1]; ts=sol.t_events[0]
    amp=float(np.max(np.abs(sol.y[0])))
    return ts,vs,amp

def strang(x,v,h,mu):
    v*=np.exp(mu*(1-x*x)*h/2.0)
    v1=v-0.5*h*x; x1=x+h*v1; v2=v1-0.5*h*x1
    v2*=np.exp(mu*(1-x1*x1)*h/2.0)
    return x1,v2

def pmap3(v0,mu,h):
    """discrete return map: start (x=0,v=v0), return v at next upward x=0 crossing"""
    x,v=0.0,v0; xp,vp=x,v
    for i in range(int(15/h)):
        x,v=strang(x,v,h,mu)
        if xp<0.0 and x>=0.0 and v>0.0:
            frac=-xp/(x-xp); return vp+frac*(v-vp)
        xp,vp=x,v
    return None

mu=0.1
ts,vs,amp=true_cycle(mu)
print("mu=0.1 TRUE cycle (Radau rtol=1e-12):  %d upward x=0 crossings"%len(ts))
print("  t at crossings:  %s"%np.array2string(ts,precision=4))
print("  v at crossings:  %s"%np.array2string(vs,precision=6))
print("  amplitude max|x| = %.8f"%amp)
v_settled=vs[-1]
print("  at last crossing: x~0, v=%.6f, F(x)~0, Liénard y=v=%.6f" %(v_settled,v_settled))
print("  Liénard-plane radius^2 = x^2+(v+F(x))^2 ~ %.6f  (circle of radius 2 => 4)"%(v_settled**2))
print("  physical energy at x=0: v^2/2 = %.6f"% (0.5*v_settled**2))
print()
print("Discrete Poincare fixed point v* (Strang, x=0 section, v>0):")
for h in [0.01,0.001,0.0005]:
    g=lambda v: pmap3(v,mu,h)-v
    vstar=brentq(g,0.5,2.5,xtol=1e-13)
    print("  h=%.4f  v*=%.8f   (Liénard y=v*=%.8f)"%(h,vstar,vstar))
print()
print("TRUE v@x=0 (settled) = %.8f ;  discrete v* (h=0.0005) -> should match."%v_settled)
