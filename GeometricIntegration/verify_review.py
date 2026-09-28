"""
verify_review.py — check the three claims in the Severe Scrutinizer review
against the actual Storer-Verlet / implicit-midpoint maps.
Run: ./.venv/Scripts/python.exe verify_review.py
"""
import sympy as sp

x, v, h = sp.symbols("x v h", real=True)
J = sp.Matrix([[0, 1], [-1, 0]])

def sympl_report(name, fx, fv):
    M = sp.Matrix([[sp.diff(fx, x), sp.diff(fx, v)],
                   [sp.diff(fv, x), sp.diff(fv, v)]])
    det = sp.simplify(M.det())
    diff = sp.simplify(M.T * J * M - J)
    zero = bool(diff.equals(sp.zeros(2, 2)))
    print("%-34s det(M) = %-14s  M^T J M - J = 0 ?  %s" % (name, str(det), zero))
    return det, zero

print("=" * 74)
print("CLAIM 1: Storer-Verlet (leapfrog) map for V = x^2/2")
print("=" * 74)
x1 = x + h*v - sp.Rational(1, 2)*h**2*x
# (a) as WRITTEN in the essay: v - hx - (h^3/4) x + (h^2/2) v
v2_essay = v - h*x - sp.Rational(1, 4)*h**3*x + sp.Rational(1, 2)*h**2*v
# (b) as COMPUTED by the script: v - hx - (h^2/2) v + (h^3/4) x
v2_script = v - h*x - sp.Rational(1, 2)*h**2*v + sp.Rational(1, 4)*h**3*x
sympl_report("(a) essay-written v2", x1, v2_essay)
sympl_report("(b) script v2 (correct)", x1, v2_script)
print("    expanded script v2 = %s" % sp.sstr(sp.expand(v2_script)))
print()

print("=" * 74)
print("CLAIM 2: Implicit midpoint (Gauss-Legendre 1-stage), V = x^2/2")
print("=" * 74)
# Midpoints (xm, vm): xm = x + (h/2) vm ; vm = v - (h/2) xm
xm = (4*x + 2*h*v)/(4 + h**2)
vm = (4*v - 2*h*x)/(4 + h**2)
# Next step (x1, v1) = 2(xm, vm) - (x, v)   [the one-step map]
X1 = ((4 - h**2)*x + 4*h*v)/(4 + h**2)
V1 = ((4 - h**2)*v - 4*h*x)/(4 + h**2)
sympl_report("midpoint map (x,v)->(xm,vm)", xm, vm)
sympl_report("next-step map (x,v)->(x1,v1)", X1, V1)
print("    next-step x1 = %s" % sp.sstr(X1))
print("    next-step v1 = %s" % sp.sstr(V1))
# check the midpoint equations are actually satisfied by (xm, vm)
print("    check xm = x + (h/2) vm :", sp.simplify(xm - (x + sp.Rational(1,2)*h*vm)) == 0)
print("    check vm = v - (h/2) xm :", sp.simplify(vm - (v - sp.Rational(1,2)*h*xm)) == 0)
# and that (X1,V1) = 2(xm,vm) - (x,v)
print("    check X1 = 2 xm - x     :", sp.simplify(X1 - (2*xm - x)) == 0)
print("    check V1 = 2 vm - v     :", sp.simplify(V1 - (2*vm - v)) == 0)
print()

print("=" * 74)
print("CLAIM 3: Radau IIA order vs. stage count")
print("=" * 74)
# s-stage Radau IIA: order p = 2s - 1, stability function = [s-1/s] Pade of exp
# s=1: order 1, R = 1/(1-z)  (backward Euler)
# s=2: order 3, R = (2 + z)/(2 - z + z^2)
z = sp.symbols("z")
R1 = 1/(1 - z)
R2 = (2 + z)/(2 - z + z**2)
print("    s=1 stage: order 2*1-1 = 1,  R(z) = 1/(1-z)   [backward Euler]")
print("    s=2 stage: order 2*2-1 = 3,  R(z) = (2+z)/(2-z+z^2)")
print("    -> 'order s' is wrong; 's-stage' (order 2s-1) is correct.")
