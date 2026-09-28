"""
Computation 5: true velocity at x=0 (v>0) on the van der Pol cycle, mu=0.1,
high accuracy (Radau, rtol=1e-12). Settles whether the discrete Poincare-map
fixed point v* -> 2.00177 converges to the TRUE section velocity or to a
crossing-detection-biased value.
Run: ./.venv/Scripts/python.exe geometric_integration_computation5.py
"""
import numpy as np
from scipy.integrate import solve_ivp

def vdp(t,y,mu):
    x,v=y; return [v, mu*(1-x*x)*v - x]

def true_section_velocity(mu, rtol=1e-12):
    # start at x=2, v=0; find upward x=0 crossings; measure v at the crossing
    def cross(t,y,*args): return y[0]
    cross.terminal=False; cross.direction=1
    # run several periods, then read v at the last few upward x=0 crossings
    t_end = 13.0   # > 2 periods (T~6.28)
    sol = solve_ivp(vdp,(0.0,t_end),[2.0,0.0],method="Radau",args=(mu,),
                    rtol=rtol,atol=rtol,events=cross,max_step=0.05)
    ye = sol.y_events[0]      # shape (2, N)
    ts = sol.t_events[0]
    vs = ye[1]                # v at each upward x=0 crossing
    amp = float(np.max(np.abs(sol.y[0])))
    return float(np.mean(vs[-5:])), amp, len(ts)

mu=0.1
vsec, amp, n = true_section_velocity(mu)
print("mu=0.1 true cycle (Radau rtol=1e-12), %d crossings in [0,13]"%n)
print("  velocity at x=0 (v>0), mean of last 5 crossings = %.8f" % vsec)
print("  max |x| (amplitude) = %.8f" % amp)
print("  (energy check: if circle of radius r, v@x=0 = r;  v^2/2 vs V(amp):",
      "%.8f vs %.8f" % (0.5*vsec**2, 0.5*amp**2))
print()
print("Discrete Poincare-map fixed point (Strang split) converged to v* = 2.0017704 (h=0.0005).")
print("  |v* - v_true| = %.3e" % abs(2.0017704 - vsec))
print("  => the discrete fixed point converges to the TRUE section velocity"
      if abs(2.0017704-vsec)<1e-4 else
      "  => WARNING: discrete v* is OFF from the true section velocity by %.3e" % abs(2.0017704-vsec))
