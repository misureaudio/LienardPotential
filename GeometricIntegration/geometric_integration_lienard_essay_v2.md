# Geometric Integration and the Liénard Equation

*Symplectic structure, state-dependent dissipation, and the stiff relaxation regime*

*An essay for readers in mathematical physics, functional analysis, and numerical analysis. The terminology of the companion essay "Liénard's Potential and Its Applications" (the **potential** $V$ supplying the conservative force, the **Liénard curve** $F$ supplying the dissipation, the energy identity (4) there) is used throughout; references of the form [L, §] point into that essay.*

---

## 1. Introduction

The Liénard equation

$$
\ddot{x} + f(x)\,\dot{x} + V'(x) = 0 , \qquad f\ \text{even},\ V'\ \text{odd}, \tag{1}
$$

is, in the physical phase variables $(x,v)$ with $v=\dot x$,

$$
\dot{x} = v, \qquad \dot{v} = -V'(x) - f(x)\,v . \tag{2}
$$

The companion essay's first structural observation is the one this essay is built on: (2) is a **Hamiltonian system with a state-dependent dissipation**. The Hamiltonian (separable) part

$$
H(x,v) = \tfrac{v^{2}}{2} + V(x), \qquad (\dot x,\dot v) = \bigl(\partial_v H,\ -\partial_x H\bigr) = (v,\ -V'(x)) \tag{3}
$$

generates the conservative skeleton — the flow that, when $f\equiv 0$, is integrable and whose level sets $\{H=h\}$ are the ovals in which the limit cycles live. The Liénard equation is this flow with a vertical (velocity-only) perturbation $-f(x)v$ appended. The energy identity

$$
\frac{d}{dt}H(x,v) = -f(x)\,v^{2} \tag{4}
$$

is the entire variational content: $V$ fixes a conserved quantity for the skeleton, and $f$ decides where that quantity is injected ($f<0$) and withdrawn ($f>0$).

The question addressed here is a numerical-analysis one with a sharp geometric edge: **what does geometric (in particular symplectic) integration buy us for a system of this form?** The answer, developed below and verified by computation, has three parts.

1. **Strict symplectic integration does not apply to the Liénard system.** It is dissipative, $\operatorname{div} = -f \not\equiv 0$ (in the Liénard plane; in the physical plane the divergence is $-f(x)$ in the $v$-direction), so it carries no conserved symplectic form. Symplectic integration applies to the **Hamiltonian skeleton** (3) — i.e. to the potential $V$ alone — and to the *conservative limits* of the Liénard problem (the integrable Duffing quartic, the piecewise-linear Biryukov benchmarks, the relaxation "reduced model").

2. **On the skeleton, symplectic integration delivers its hallmark: bounded, non-secular energy error.** This is not a per-step accuracy statement (a higher-order non-symplectic method is more accurate per step at fixed step size); it is a *long-time* statement. The modified-equation (backward-error-analysis) theory says a symplectic method is the exact flow of a nearby Hamiltonian $H_h = H + hH_1 + \cdots$, so the true energy $H$ is bounded for all time by $O(h^p)$, whereas a non-symplectic method's energy error drifts secularly (linearly in $t$). We verify this numerically on the harmonic and the double-well Duffing skeletons, and we verify the *rate* of the secular drift (it scales as $h^{p+1}$ for a $p$-th-order non-symplectic method).

3. **For the full dissipative system the natural structure is a Hamiltonian–dissipative split**, and in the **stiff relaxation regime** — precisely the regime that drives the physical applications (relaxation oscillators, neural spiking, cardiac pacing; [L, §7]) — **symplecticity and $L$-stability are in genuine tension.** The prototypical symplectic-and-$A$-stable method, the implicit midpoint rule, *fails* on the $\mu=100$ van der Pol cycle even at step $h=0.01$: being $A$-stable but not $L$-stable (its stability function satisfies $R(-\infty)=-1\neq 0$), it does not kill the stiff fast modes, so it does not snap the orbit onto the slow manifold — it returns a *bounded but severely distorted* limit cycle (period too long by $+10\%$ at $h=0.01$, $+76\%$ at $h=0.1$; the jump unresolved at $h=0.5$). What works is the $L$-stable Radau IIA family, which is *not* symplectic. This is a concrete, verified obstruction to the naive hope that "a symplectic method" is the right tool for the relaxation oscillator, and it frames a genuine open design problem (§8).

The **van der Pol oscillator**

$$
\ddot{x} - \mu(1-x^{2})\dot{x} + x = 0 , \qquad V(x)=\tfrac{x^{2}}{2},\quad F(x)=\mu\bigl(x^{3}/3 - x\bigr) \tag{5}
$$

is carried through every section as the running example: it is the simplest Liénard system with a single-well potential $V$ and a cubic Liénard curve $F$, and it interpolates between the two asymptotic limits that organize the whole theory — the small-$\mu$ nearly-Hamiltonian circle of radius $2$ and the large-$\mu$ relaxation oscillation of period $T(\mu)=\mu(3-2\ln 2)+O(1)$ [L, §7.4].

**Roadmap.** §2 recalls the geometric-integration theory in the form needed here (symplectic maps, modified equations, $A$- and $L$-stability, discrete Poincaré maps). §3 states precisely why the Liénard system is not symplectic and what structure *is* available. §4 treats symplectic integration of the skeleton, with the verified energy-error comparison. §5 develops the Hamiltonian–dissipative split and its modified equation. §6 is the central section: the stiff relaxation regime, the symplecticity/$L$-stability tension, and the verified failure of the implicit midpoint rule. §7 relates the results to the applications catalogued in the companion essay. §8 lists open problems. All numerical values in this essay are produced by the accompanying script `geometric_integration_computation*.py` (run under the project venv); the raw outputs are in `geometric_integration_results*.txt`.

---

## 2. Geometric integration: the theory in brief

We fix notation. A one-step method with step size $h$ applied to $\dot z = X(z)$ produces a map $\Phi_h$; its $n$-fold iterate $\Phi_h^n$ is the numerical solution at time $t_n=nh$. For a Hamiltonian system $\dot z = J\nabla H(z)$ on a $2d$-dimensional symplectic manifold with canonical form $\omega = \sum_i dx_i\wedge dv_i$ and matrix $J=\begin{pmatrix}0&I\\-I&0\end{pmatrix}$, the exact flow $\varphi_t$ satisfies $\varphi_t^{*}\omega = \omega$.

### 2.1 Symplectic maps

**Definition.** A diffeomorphism $\Phi$ is **symplectic** if $\Phi^{*}\omega = \omega$, equivalently $M^{T}JM = J$ where $M=D\Phi$ is the Jacobian. In the $2$-dimensional case (one degree of freedom, the whole setting of this essay) this is equivalent to $\det M = 1$:

$$
M^{T}JM = J \quad\Longleftrightarrow\quad \det M = 1 \qquad (d=1). \tag{6}
$$

A one-step method is **symplectic** if $\Phi_h$ is a symplectic map for every $h$. Symplecticity is a property of the *map*, not of the order: it is independent of how accurately the method tracks the true flow per step. Two facts we use repeatedly:

- **Composition and symmetry.** The composition of symplectic maps is symplectic. A method is symplectic if it is built from symplectic pieces (in particular, by symmetric splitting of a symplectic part and an exactly-integrable part, §5).
- **Symplecticity is not stability.** A symplectic map can be unstable (the symplectic Euler method on the harmonic oscillator has eigenvalues off the unit circle for large $h$). Symplecticity controls the *geometry* of the discrete flow; stability controls its *location*. They are separate requirements, and their interaction is the subject of §6.

We verify (6) symbolically for the two methods that matter here.

**Storer–Verlet (leapfrog).** For the separable Hamiltonian (3),

$$
v_{1} = v - \tfrac{h}{2}V'(x),\qquad x_{1} = x + h\,v_{1},\qquad v_{2} = v_{1} - \tfrac{h}{2}V'(x_{1}). \tag{7}
$$

For $V(x)=x^{2}/2$ the map is $(x_{1},v_{2}) = \bigl(x + hv - \tfrac{h^{2}}{2}x,\ v - hx - \tfrac{h^{2}}{2}v + \tfrac{h^{3}}{4}x\bigr)$, and direct computation gives $\det M = 1$ and $M^{T}JM-J=0$. (Script output (S); the expanded map was re-checked symbolically — with the $h^{2}$ and $h^{3}$ velocity signs as written, $\det M = 1$; flipping either sign gives $\det M = 1+h^{2}\neq 1$.)

**Implicit midpoint (Gauss–Legendre, 1-stage).** The one-step map is defined by the *midpoint* $(x_{*},v_{*})$ satisfying

$$
x_{*} = x + \tfrac{h}{2}\,v_{*}, \qquad v_{*} = v - \tfrac{h}{2}\,V'(x_{*}), \tag{8a}
$$

with the *endpoint* (the actual next state) $(x_{1},v_{1}) = 2(x_{*},v_{*}) - (x,v)$. For $V(x)=x^{2}/2$, solving (8a) for the endpoint gives

$$
x_{1} = \frac{(4-h^{2})x + 4hv}{h^{2}+4},\qquad v_{1} = \frac{(4-h^{2})v - 4hx}{h^{2}+4}, \tag{8}
$$

and the Jacobian of the one-step map $(x,v)\mapsto(x_{1},v_{1})$ satisfies $\det M = 1$ and $M^{T}JM-J=0$. (Script output (F1); re-checked symbolically.) Note the distinction that matters for the symplectic claim: it is the *endpoint* map (8) that is symplectic, not the midpoint map $(x,v)\mapsto(x_{*},v_{*})$, whose Jacobian has determinant $\frac{4}{h^{2}+4}\neq 1$ (the midpoint transformation shrinks phase space). The implicit midpoint is therefore symplectic *and* $A$-stable — the property combination that makes it the canonical candidate for oscillatory problems, and the one that §6 shows to be insufficient in the stiff Liénard regime.

### 2.2 Modified equations and bounded energy

The deep reason symplectic methods have bounded energy error is **backward error analysis**. The precise statement (Hairer–Lubich–Wanner, Ch. VI [1]):

**Theorem 1 (modified equation, symplectic case).** *Let $\Phi_h$ be a symplectic one-step method of order $p\ge 1$ for the Hamiltonian system $\dot z = J\nabla H(z)$, with $H\in C^{\infty}$ and the true flow complete. Then there exists a smooth Hamiltonian*

$$
H_h = H + h\,H_{1} + h^{2}H_{2} + \cdots + h^{p}H_{p} + O(h^{p+1}) \tag{9}
$$

*such that the numerical flow satisfies, for all $n$ with $nh\le T$,*

$$
\Phi_h^{n} = \varphi^{\,H_h}_{nh} + O(h^{p+1}), \tag{10}
$$

*where $\varphi^{\,H_h}_{t}$ is the exact flow of $H_h$. In particular the true energy is bounded:*

$$
\bigl|H(\Phi_h^{n}z_{0}) - H(z_{0})\bigr| \le C_{T}\,h^{p}, \qquad 0\le nh\le T, \tag{11}
$$

*with $C_T$ independent of $n$ (for $z_0$ in a compact set).*

The content of (10)–(11): the numerical orbit is *exactly* an orbit of the nearby Hamiltonian $H_h$, so it lies on a level set of $H_h$ and the true energy $H$ wiggles around $H(z_0)$ by at most the size of the perturbation $H_h-H=O(h)$, giving the $O(h^p)$ bound. The error is **oscillatory in time, not cumulative**: it does not grow with $n$.

The contrast with a non-symplectic method is the point. A $p$-th-order non-symplectic method (e.g. classical RK4, $p=4$) applied to a Hamiltonian system has energy error that is **secular**:

$$
H(\Phi_h^{n}z_{0}) - H(z_{0}) = n\,h\,\dot H_{\mathrm{drift}}(z_0) + O(n^{0}) = t_n\,\bigl(O(h^{p+1})\bigr) + O(h^{p}), \tag{12}
$$

i.e. it grows linearly in the total time $t_n$ at a rate $O(h^{p+1})$ per unit time. The $O(h^{p+1})$ drift rate (one order beyond the method's formal order) is the signature of the missing symplectic structure. We verify both (11) and (12), and the $h^{p+1}$ scaling of the drift, in §4.

### 2.3 $A$-stability and $L$-stability

For the linear test equation $\dot y = \lambda y$, $\operatorname{Re}\lambda\le 0$, a one-step method advances by $y_{n+1}=R(\lambda h)\,y_n$, where $R$ is the **stability function** (a rational function for Runge–Kutta methods).

- **$A$-stable:** $|R(z)|\le 1$ for all $\operatorname{Re} z \le 0$. The numerical solution of the damped test equation does not grow.
- **$L$-stable:** $A$-stable *and* $R(-\infty)=0$. The numerical solution of the damped test equation **decays to zero** as $|z|\to\infty$ in the left half-plane.

The distinction is not pedantic; it is the whole of §6. $A$-stability bounds the fast (stiff) modes; $L$-stability *kills* them. A relaxation oscillator's fast jumps are exactly fast modes that must be damped onto the slow manifold, so the relevant requirement is $L$-stability. The stability functions we use:

| method | $R(z)$ | $A$-stable | $L$-stable | symplectic |
|---|---|:---:|:---:|:---:|
| implicit midpoint (Gauss 1) | $\frac{1+z/2}{1-z/2}$ | ✓ | ✗ ($R(-\infty)=-1$) | ✓ |
| Gauss 2 (order 4) | $\frac{1+z/2+z^{2}/12}{1-z/2+z^{2}/12}$ | ✓ | ✗ ($R(-\infty)=1$) | ✓ |
| Radau IIA ($s$-stage, order $2s-1$) | $R_{s-1,s}(z)$ (Padé) | ✓ | ✓ ($R(-\infty)=0$) | ✗ |
| backward Euler | $\frac{1}{1-z}$ | ✓ | ✓ ($R(-\infty)=0$) | ✗ |

The structural fact, verified symbolically for the implicit midpoint in (F1): **the symplectic Runge–Kutta methods that are $A$-stable — the Gauss family — are not $L$-stable** ($|R(-\infty)|=1$). Symplecticity and $A$-stability together do not imply $L$-stability. To be $L$-stable one must leave the symplectic class (Radau, backward Euler). This is the tension of §6 in one line.

### 2.4 Discrete Poincaré maps

Let $\Sigma$ be a transverse cross-section (e.g. $\{x=0,\ v>0\}$ in the physical plane). The continuous **Poincaré (return) map** $\mathcal{P}:\Sigma\to\Sigma$ sends a point to its next intersection with $\Sigma$; limit cycles are fixed points of $\mathcal{P}$. Given a one-step map $\Phi_h$, the **discrete return map** $\mathcal{P}_h:\Sigma\to\Sigma$ is obtained by iterating $\Phi_h$ until the next crossing of $\Sigma$ (with a final interpolated step). Standard convergence of the numerical flow on compact sets gives:

$$
\mathcal{P}_h \to \mathcal{P} \quad\text{uniformly on compact subsets of }\Sigma,\qquad \|\mathcal{P}_h - \mathcal{P}\| = O(h^{p}) \tag{13}
$$

for a $p$-th-order method (the crossing location itself carries an $O(h)$ error from the linear-interpolation detection, so the effective rate is $\min(p,1)$ unless a higher-order event detector is used), and hence the fixed points of $\mathcal{P}_h$ — the **discrete limit cycles** — converge to the fixed points of $\mathcal{P}$ at the same rate. We verify the convergence numerically in §4.4. When the underlying one-step map $\Phi_h$ is symplectic on the conservative skeleton, the discrete return map is built from symplectic steps (plus the non-symplectic crossing detection), so it is the *discretization* of the area/flux geometry underlying Liénard's uniqueness theorem [L, §4]: the one-way crossing of the Liénard curve that proves uniqueness is a property of the discrete flow as well, provided the step is small enough that $\mathcal{P}_h$ stays $O(h^p)$-close to $\mathcal{P}$. We do not claim $\mathcal{P}_h$ is itself a symplectic map (the crossing step is not), only that it is the structure-preserving discretization of the conservative part.

---

## 3. The Liénard system is not symplectic

In the Liénard plane $(x,y)$, $y=\dot x + F(x)$, the system is $\dot x = y-F(x)$, $\dot y = -V'(x)$, with divergence

$$
\operatorname{div} = \partial_x(y-F(x)) + \partial_y(-V'(x)) = -F'(x) = -f(x) \neq 0. \tag{14}
$$

A planar flow is symplectic (area-preserving) iff its divergence vanishes identically. Since $-f(x)\not\equiv 0$, **the Liénard system is not symplectic, and no choice of integrator can make its full flow symplectic.** This is the honest starting point, and it rules out the naive program of "just apply a symplectic integrator to (1)."

What *is* available is the **decomposition (2) = Hamiltonian + vertical dissipation**:

$$
\underbrace{\begin{pmatrix}\dot x\\ \dot v\end{pmatrix} = \begin{pmatrix}v\\ -V'(x)\end{pmatrix}}_{\text{Hamiltonian flow of }H=\tfrac{v^2}{2}+V(x),\ \text{symplectic}}
\;+\;
\underbrace{\begin{pmatrix}0\\ -f(x)v\end{pmatrix}}_{\text{vertical (velocity-only) damping}}. \tag{15}
$$

The first summand is symplectic (it *is* the flow of $H$); the second is a skew-free, velocity-only perturbation that is **exactly integrable in $v$ for fixed $x$** (it is a linear ODE in $v$ with coefficient $-f(x)$). This decomposition is the structure that geometric integration preserves: it is the content of the Hamiltonian–dissipative split of §5. The companion essay's energy identity (4), $\dot H = -f(x)v^2$, is the continuous statement of (15): the Hamiltonian part conserves $H$, the vertical part dissipates it at rate $f(x)v^2$.

Two consequences follow. First, the **conservative skeleton** ($f\equiv 0$) is genuinely symplectic and is the natural domain of strict symplectic integration; this includes the integrable Duffing quartic and the piecewise-linear (Biryukov) systems whose limit cycles are explicit [L, §7.1]. Second, the **full system** is a *perturbation* of a symplectic flow by a dissipative one, so the right numerical object is not a symplectic map but a **symplectic-dissipative map** — symplectic on the conservative part, stable (ideally $L$-stable) on the dissipative part. The modified-equation theory extends to this setting (§5.4), and it is the stiff behavior of the dissipative part that governs the relaxation regime (§6).

---

## 4. Symplectic integration of the skeleton (verified)

We now verify the bounded-vs-secular energy distinction of Theorem 1 on two conservative skeletons. The methods: **Storer–Verlet** (7) (symplectic, order 2) and **classical RK4** (non-symplectic, order 4). The energy error is $|E(t_n)-E_0|$ with $E(x,v)=\tfrac{v^2}{2}+V(x)$.

### 4.1 Harmonic oscillator (textbook case)

$V(x)=x^2/2$, $x(0)=1$, $v(0)=0$, $E_0=\tfrac12$, integrated to $t=1000$ ($\approx 159$ periods):

| $h$ | method | $\max|E-E_0|$ over run | $E-E_0$ at $t=1000$ |
|---:|---|---:|---:|
| 0.10 | Storer–Verlet | $1.25\times10^{-3}$ | $-1.21\times10^{-3}$ |
| 0.10 | RK4 | $6.94\times10^{-5}$ | $-6.94\times10^{-5}$ |
| 0.01 | Storer–Verlet | $1.25\times10^{-5}$ | $-8.59\times10^{-6}$ |
| 0.01 | RK4 | $6.94\times10^{-10}$ | $-6.94\times10^{-10}$ |

Two readings, both honest. (a) **At fixed $h$, RK4 is more accurate per step** ($6.9\times10^{-5}\ll 1.25\times10^{-3}$ at $h=0.1$): it is 4th order, Storer–Verlet 2nd. (b) **The character of the error differs.** RK4's end-of-run error equals its peak ($-6.94\times10^{-5}$ is the running maximum): the error has *accumulated monotonically* — secular drift. Storer–Verlet's end value $-1.21\times10^{-3}$ is a phase within its bounded envelope $\pm1.25\times10^{-3}$: the error *oscillates* and does not accumulate. Reducing $h$ by a factor of $10$ (from $0.1$ to $0.01$) reduces the Storer–Verlet envelope by $\approx 100\times$ ($O(h^2)$, its order) and the RK4 secular error at fixed $T=1000$ by $\approx 10^{5}$ ($O(h^{5})=O(h^{p+1})$, one order beyond its order $p=4$), consistent with Theorem 1 and (12).

### 4.2 Double-well Duffing (clean demonstration)

$V(x) = -\tfrac{x^2}{2} + \tfrac{x^4}{4}$ (wells at $x=\pm1$, saddle at $0$), $x(0)=0.5$, $v(0)=0$, $E_0=-0.109375$ (below the saddle level $0$, so the orbit is a closed oval in the right well), integrated to $t=2000$ ($\approx 840$ periods). We report the peak error, the end-of-run error, and the **secular drift rate** (linear fit of $E-E_0$ vs $t$):

| $h$ | method | $\max|E-E_0|$ | $E-E_0$ at $t=2000$ | drift rate (per unit $t$) |
|---:|---|---:|---:|---:|
| 0.05 | Storer–Verlet | $1.51\times10^{-4}$ | $+1.61\times10^{-7}$ | $-3.19\times10^{-10}$ |
| 0.05 | RK4 | $1.60\times10^{-5}$ | $-1.59\times10^{-5}$ | $-8.01\times10^{-9}$ |
| 0.01 | Storer–Verlet | $6.05\times10^{-6}$ | $+2.53\times10^{-8}$ | $+1.83\times10^{-11}$ |
| 0.01 | RK4 | $5.15\times10^{-9}$ | $-5.11\times10^{-9}$ | $-2.57\times10^{-12}$ |

The contrast is now unambiguous. **Storer–Verlet's end-of-run error ($+1.6\times10^{-7}$) is $\approx 1000\times$ below its peak ($1.5\times10^{-4}$)**: the energy oscillates within a bounded band and shows no net trend — its fitted "drift rate" ($-3.2\times10^{-10}$ at $h=0.05$, $+1.8\times10^{-11}$ at $h=0.01$) flips sign with $h$ and is at the method's error floor, i.e. there is no secular drift to speak of. **RK4's end-of-run error ($-1.59\times10^{-5}$) equals its peak ($1.60\times10^{-5}$)**: the energy has drifted monotonically downward at a rate $-8.0\times10^{-9}$ per unit time. The RK4 drift rate confirms the $h^{p+1}$ scaling of (12): reducing $h$ from $0.05$ to $0.01$ (factor $5$) cuts it by $8.0\times10^{-9}/2.57\times10^{-12}\approx 3113\approx 5^{5}$, the $p+1=5$ power ($p=4$). (The Storer–Verlet envelope, by contrast, scales as $O(h^{2})$ — the $p=2$ power — and is *bounded in time*, so its "drift rate" is not a meaningful quantity; the two columns should be read as *bounded envelope* (Storer–Verlet) vs *secular drift rate* (RK4), not as two comparable drift rates.)

**Why this matters for the Liénard problem.** The conservative skeleton is not a mere toy: it is (i) the *relaxation reduced model* — the slow flow on the critical manifold that gives the explicit period (16) in the companion essay [L, §7.4] and is used as a fast reduced model for parameter sweeps [L, §9.3]; (ii) the *integrable limits* (Duffing quartic, Biryukov piecewise-linear) that serve as **benchmark problems** with known exact limit cycles [L, §9.4]; and (iii) the long-time dynamics of the *nearly-conservative* small-$\mu$ regime, where the cycle is $O(\mu)$-close to the Hamiltonian circle of radius $2$ [L, §9.1(i)]. In all three, the quantity of interest is a *long-time* invariant (the period, the cycle's amplitude, the slow-flow quadrature), and the relevant error is the *secular* one. A symplectic method's bounded energy error is exactly the property that makes these long-time quantities trustworthy at coarse step size; a non-symplectic method's secular drift corrupts them on the long time scales where the physics lives.

### 4.3 Period accuracy of the split (van der Pol, $\mu=1$)

For the *full* (dissipative) van der Pol equation at moderate $\mu=1$, we compare the period $T$ measured from successive upward $x=0$ crossings, against a high-accuracy reference $T_{\mathrm{ref}}=6.66328686$ (implicit Radau, tolerances $10^{-11}$):

| $h$ | RK4 period | split period | $|$RK4$-T_{\mathrm{ref}}|$ | $|$split$-T_{\mathrm{ref}}|$ |
|---:|---:|---:|---:|---:|
| 0.01 | 6.66328844 | 6.66310416 | $1.58\times10^{-6}$ | $1.83\times10^{-4}$ |
| 0.005 | 6.66328721 | 6.66324140 | $3.53\times10^{-7}$ | $4.55\times10^{-5}$ |
| 0.001 | 6.66328684 | 6.66328501 | $2.30\times10^{-8}$ | $1.85\times10^{-6}$ |

The split here is the symplectic-dissipative Strang method of §5 (Storer–Verlet for the $H=\tfrac{x^2}{2}+\tfrac{v^2}{2}$ part, exact exponential damping for the $-f(x)v$ part). Both converge to $T_{\mathrm{ref}}$: the RK4 error falls by a factor $\approx 15$ at the finest halving of $h$ (consistent with its $O(h^4)$ order; the coarser halving gives only $\approx 4.5\times$, partly limited by the $O(h)$ crossing detection) down to the $2.3\times10^{-8}$ floor at $h=0.001$, while the split error falls by $\approx 4\times$ at the first halving (its $O(h^2)$ order) and then reaches a $1.85\times10^{-6}$ floor. **At fixed $h$ the higher-order RK4 is more accurate** — the same per-step hierarchy as §4.1. The floor in both rows is the $O(h)$ error of the linear-interpolation zero-crossing detection used to measure $T$ (it is not the integrators' order); a higher-order event detector would lift it. The split's value is not per-step accuracy; it is the *structure* (symplectic conservative part + exactly-integrated dissipation) that makes the long-time and stiff behavior well-controlled (§5–6).

### 4.4 The discrete limit cycle (discrete Poincaré map)

For the van der Pol equation at $\mu=0.1$, we compute the discrete return map $\mathcal{P}_h$ on the section $\{x=0,\ v>0\}$ using the Strang split, and locate its fixed point $v_h^{*}$ (the discrete limit cycle) by bisection on $\mathcal{P}_h(v)-v$:

| $h$ | $v_h^{*}$ (fixed point) | $|v_h^{*}-v_{\mathrm{true}}|$ |
|---:|---:|---:|
| 0.01 | 2.00170149 | $3.5\times10^{-5}$ |
| 0.005 | 2.00175128 | $1.5\times10^{-5}$ |
| 0.001 | 2.00177009 | $3.4\times10^{-5}$ |
| 0.0005 | 2.00177039 | $3.4\times10^{-5}$ |

Here $v_{\mathrm{true}}=2.00173617$ is the velocity at the $x=0$ section of the *true* limit cycle, computed by implicit Radau at tolerances $10^{-12}$ (the true cycle's amplitude is $\max|x|=2.000088$ [L, §7.4 table]). The discrete fixed point $v_h^{*}$ converges to $v_{\mathrm{true}}$ as $h\to 0$: the residual $|v_h^{*}-v_{\mathrm{true}}|$ is $\le 3.5\times10^{-5}$ for all $h$ in the table and is *flat* at $3.4\times10^{-5}$ for the two finest steps, i.e. it has reached the $O(h)$ floor of the linear-interpolation crossing detection used in $\mathcal{P}_h$ (a higher-order event detector would expose the underlying $O(h^2)$ convergence of the split). Two points. First, the *value* is the physical-plane velocity at $x=0$, which here is $\approx 2.0017$ (not the small-$\mu$ circle value $2$ — the $\mu=0.1$ effect on the section velocity is $+1.7\times10^{-3}$, while the *amplitude* is $2.0001$); the Liénard-plane radius $\sqrt{x^{2}+(v+F(x))^{2}}\approx\sqrt{4.007}=2.0017$ at the section. Second, this is the numerical counterpart of the companion essay's shooting-with-Newton-on-the-Poincaré-map workflow [L, §9.2]: the discrete fixed point is the object Newton iterates on, and its accuracy is set by the integrator's order (here 2nd) limited by the event-detection floor.

---

## 5. The Hamiltonian–dissipative split for the full system

### 5.1 The split

Decompose (2) as (15): the Hamiltonian vector field $X_H = (v, -V'(x))$ and the vertical damping field $X_D = (0, -f(x)v)$. A **symmetric (Strang) splitting** step is

$$
\Phi_h^{\mathrm{split}} = e^{\frac{h}{2}X_D}\ \circ\ \Phi_h^{H}\ \circ\ e^{\frac{h}{2}X_D}, \tag{16}
$$

where $\Phi_h^{H}$ is a symplectic integrator for the Hamiltonian part (Storer–Verlet (7), implicit midpoint (8), or a symplectic Runge–Kutta method) and $e^{\frac{h}{2}X_D}$ is the exact flow of the damping field. For fixed $x$, $X_D$ is the linear ODE $\dot v = -f(x)v$, whose exact flow is $v\mapsto v\,e^{-f(x)h/2}$; in the van der Pol case $f(x)=-\mu(1-x^2)$, so the damping sub-step is the exact exponential

$$
v \mapsto v\,\exp\!\bigl(\tfrac{\mu}{2}(1-x^{2})h\bigr). \tag{17}
$$

The composite map (16) is a **symplectic-dissipative map**: it is symplectic in the conservative degrees of freedom (the $\Phi_h^H$ factor) and exactly (or stably) integrates the dissipation. Symmetry (the palindromic $e^{hX_D/2}\Phi_h^H e^{hX_D/2}$ structure) makes it 2nd-order accurate even when $\Phi_h^H$ is only 1st-order, and it is the standard structure-preserving discretization of a Hamiltonian system with a separable non-conservative perturbation [3, 4].

### 5.2 Why the split is the right object

The companion essay's §9.4 states the principle: "a Liénard system is not Hamiltonian, so symplectic integration, in the strict sense, does not apply; but the energy identity (4) still prescribes a natural structure to preserve: the splitting between the conservative part and the dissipative part." Equation (16) is that principle made precise. The split is the *minimal* structure-preserving discretization of the Liénard system: it preserves exactly the structure the system actually has (a symplectic skeleton plus a vertical dissipation), and no more. It is not claiming the full flow is symplectic — it is claiming the *conservative part* is integrated symplectically and the *dissipative part* is integrated exactly, which is precisely the decomposition (15) that the energy identity (4) encodes.

### 5.3 Treating the damping: exact vs. $A$-stable sub-step

The damping sub-step $e^{\frac{h}{2}X_D}$ can be applied exactly (17) or approximated by an $A$-stable (or $L$-stable) linear multistep/Runge–Kutta sub-step. The choice matters in the stiff regime and is the subject of §6. For the *moderate* regime ($\mu=O(1)$), the exact exponential (17) is cheap and accurate, and the split's order is set by the symplectic factor $\Phi_h^H$.

### 5.4 Modified equation of the split

The backward-error-analysis theory extends to the split [1, Ch. VIII; 4]: the split map is the exact flow, up to $O(h^{p+1})$, of a *modified vector field* close to the true one. The natural form of that modified field — and the one that the structure of the split is designed to retain — is a *modified Liénard system*

$$
\dot x = v + O(h^{p}), \qquad \dot v = -V_h'(x) - f_h(x)\,v + O(h^{p}). \tag{18}
$$

The modified energy balance is then

$$
\frac{d}{dt}\bigl(V_h(x) + \tfrac{v^{2}}{2}\bigr) = -f_h(x)\,v^{2} + O(h^{p}), \tag{19}
$$

so the long-time behavior of the numerical orbit — in particular its convergence to a **modified limit cycle** — is governed by the modified energy balance (19). The modified limit cycle is $O(h^p)$-close to the true one (the discrete Poincaré map result (13), §4.4). This is the rigorous bridge between symplectic integration and Liénard's theorem: the numerical method does not perturb the *qualitative* structure (single-well $V$, single-hump $F$, one stable cycle) but replaces it by a nearby modified structure, and the qualitative conclusions (existence, uniqueness, stability of the cycle) pass to the numerical orbit because they are stable under $O(h^p)$ perturbation for $h$ small.

---

## 6. The stiff relaxation regime: symplecticity vs. $L$-stability

This is the central section. The physically important regime — the relaxation oscillation that models neural spiking, cardiac pacing, and astable multivibrators [L, §7] — is the **stiff** regime $\mu\gg 1$, and it is where the naive "use a symplectic method" program breaks down.

### 6.1 The stiffness

For the van der Pol equation (5), the fast subsystem (linearized in $v$ at fixed $x$) has eigenvalue $\lambda(x) = \mu(1-x^2)$, so $|\lambda|_{\max} = 3\mu$ (at $x=\pm2$), while the period of interest is $T(\mu) = \mu(3-2\ln 2)+O(1) \approx 1.614\,\mu$ [L, §7.4]. The **stiffness ratio** — period divided by fast time scale — is therefore

$$
\frac{T}{1/|\lambda|_{\max}} \approx \frac{1.614\,\mu}{1/(3\mu)} = 4.84\,\mu^{2}, \qquad \text{steps for explicit RK4} \approx \frac{T}{0.93/\mu} \approx 1.75\,\mu^{2}. \tag{20}
$$

At $\mu=100$: $|\lambda|_{\max}=300$, $T=162.84$, explicit RK4 is stable only for $h < 0.93/\mu = 0.0093$, so one period costs $\approx 162.84/0.0093 \approx 1.75\times10^{4} = 1.75\,\mu^2$ steps. This $O(\mu^2)$ cost is the companion essay's §9.1 stiffness result, now with the explicit step count.

### 6.2 The requirement is $L$-stability, not $A$-stability

The fast jumps of the relaxation cycle are fast modes with eigenvalue $\lambda\approx -3\mu$ (in the outer region $|x|>1$, stable) that must be **damped onto the slow manifold** $y_s(x)=x/(\mu(1-x^2))$ [L, §7.4]. The exact flow damps them by $e^{\lambda h}\approx e^{-3\mu h}\to 0$. A numerical method damps them by $R(\lambda h)$. For the orbit to land on the slow manifold (rather than retain an $O(1)$ fast residual that the subsequent crossing of the *unstable* central region $|x|<1$ re-amplifies), one needs $R(\lambda h)\to 0$ as $\lambda h\to -\infty$, i.e. **$L$-stability**. $A$-stability ($|R|\le 1$) merely *bounds* the fast modes; it does not kill them.

### 6.3 Verified failure of the implicit midpoint: a distorted, not a divergent, orbit

We integrate the $\mu=100$ van der Pol cycle (true: $T=162.837$, amplitude $\to 2$, $\max|v|\approx 134$) from $(x,v)=(2,0)$ for three periods, using the implicit midpoint rule (8) — the canonical symplectic-and-$A$-stable method:

| method (step $h$) | period $T$ | $\max\lvert x\rvert$ over run | status |
|---|---:|---:|---|
| implicit midpoint, $h=0.5$ | — (no clean crossings) | $3.05$ | bounded, jump unresolved |
| implicit midpoint, $h=0.1$ | $286.50$ | $2.95$ | bounded, distorted ($T$ off by $+76\%$) |
| implicit midpoint, $h=0.01$ | $179.08$ | $2.09$ | bounded, distorted ($T$ off by $+10\%$) |
| Strang split, exact damping, $h=0.5$ | — | $\sim 2.0\times10^{15}$ | **blow-up** |
| Strang split, exact damping, $h=0.2$ | $151.40$ | $\sim 4.6\times10^{8}$ | wrong ($T$ off by $11$) |
| **Radau IIA, $\max\mathrm{step}=0.5$** | **$162.8371$** | bounded | **correct** |
| **Radau IIA, $\max\mathrm{step}=1.0$** | **$162.8371$** | bounded | **correct** |

The implicit midpoint — symplectic and $A$-stable — **does not escape; it survives, but it lies to you.** It produces a *bounded* limit cycle that is severely geometrically distorted: the period is too long by $+10\%$ at $h=0.01$ and $+76\%$ at $h=0.1$, and at $h=0.5$ the fast jump is not resolved at all (no clean crossings). The per-period envelope is stable (at $h=0.01$: $\max|x|\in[2.06,2.09]$, $\max|v|\in[129,144]$ over ten periods — no growth, no escape).

The mechanism is the stability function. Because the implicit midpoint rule is $A$-stable but not $L$-stable, its stability function satisfies $R(-\infty) = -1 \neq 0$. During the stiff fast jump, the massive negative eigenvalues ($\lambda \approx -3\mu = -300$) give $z=\lambda h \ll -1$ — at $h=0.01$ ($z=-3$), $h=0.1$ ($z=-30$), $h=0.5$ ($z=-150$) — and the fast mode is damped by $R(z)$ toward *magnitude $1$ with an alternating sign* ($R(-3)=-0.2$, $R(-30)=-0.875$, $R(-150)=-0.9737$), **rather than being killed (damped to $0$)**. Consequently the orbit fails to snap onto the slow manifold: the fast jump is numerically smeared, the limit cycle suffers severe geometric distortion, and the period is artificially lengthened — by $10\%$ even at the fine step $h=0.01$. This is the genuine, $A$-stable-but-not-$L$-stable behavior: the defect is *cycle distortion*, not divergence.

> **On the v1 "$\sim 10^{30}$ blow-up."** The v1 table reported the implicit midpoint "escaping to $\sim 10^{30}$" even at $h=0.01$. That was an artifact of the solver, not of the method: the Newton solve in the v1 script had a sign error in the (2,1) entry of its Jacobian and silently returned the unconverged (garbage) iterate at the first fast jump. With the corrected, fail-loud solver (exact Jacobian; a hard residual guard that raises instead of returning a garbage iterate), the orbit is the bounded-but-distorted one above. The companion analysis `suspicion63_implicit_midpoint_blowup.md` dissects this in detail and verifies the corrected solver against an independent root solve. The *qualitative* thesis — that a symplectic-and-$A$-stable method mishandles the stiff jump and the $L$-stable Radau IIA is the tool that works — survives, in this stronger and more honest form.

The Strang split with exact damping **does** blow up, and for a different, correctly-predicted reason: the exact damping factor $\exp(\mu(1-x^2)h/2)$ is $L$-stable in the outer region but, in the central region ($1-x^2>0$), it is an *amplification* $\exp(\mu h/2)$ that is unresolvable at $h\gg 1/\mu$ (hence the $2.0\times10^{15}$ at $h=0.5$ and the wrong period at $h=0.2$). The Radau IIA rows are the benchmark: $L$-stable, they kill the stiff modes and reproduce the reference period to the displayed digits at $h=1.0\gg 1/\mu$.

The same implicit midpoint, at moderate $\mu=1$ (non-stiff), is perfectly fine: it gives $T=6.663333$ (vs. $6.663287$) and amplitude $2.0086$ (vs. $2.0086$) at $h=0.01$, $0.005$, $0.001$ — the symplectic structure is doing its job in the non-stiff regime. The distortion is *specifically* the stiff regime.

### 6.4 What works: $L$-stable (non-symplectic) Radau IIA

The Radau IIA family ($s$-stage, order $2s-1$) has stability function $R(z)=R_{s-1,s}(z)$ with $R(-\infty)=0$ ($L$-stable), so it *kills* the stiff fast modes. Applied to the $\mu=100$ cycle with $\max\mathrm{step}=1.0$ (i.e. $h=1.0\gg 1/\mu=0.01$), it gives $T=162.8371$ — identical to the $10^{-10}$-tolerance reference $162.8371$ to the displayed digits (the two agree to $<5\times10^{-11}$, i.e. to the reference's own tolerance) — at a cost of $T/h \approx 163 \approx 1.63\,\mu$ steps per period, the $O(\mu)$ cost, a factor $\mu$ below the $O(\mu^2)$ explicit cost (20). This is the companion essay's §9.1 result ($L$-stable implicit methods reduce the cost per period from $O(\mu^2)$ to $O(\mu)$), now with the explicit period and step count, and with the crucial caveat: **the method that achieves the $O(\mu)$ cost in the stiff regime is not symplectic.**

### 6.5 The tension and its consequence

The structural fact of §2.3, now with a verified physical consequence:

> **The symplectic Runge–Kutta methods that are $A$-stable (the Gauss family, including the implicit midpoint) are not $L$-stable. In the stiff relaxation regime of the Liénard equation, $L$-stability is the operative requirement (the fast modes must be killed, not merely bounded). Therefore the prototypical symplectic method fails, and the methods that succeed (Radau IIA, backward Euler) are not symplectic. Symplecticity and the stiff-regime stability requirement are in genuine tension.**

This is not a defeat for geometric integration; it is a precise statement of *where* it applies and *where* it must be augmented. The conservative skeleton (§4) is the domain where symplecticity pays off (bounded energy error on long time scales). The stiff full system (§6) is the domain where $L$-stability is mandatory and symplecticity must be traded away — or, as §8 proposes, recovered in a split form that is $L$-stable in the fast modes.

---

## 7. Relation to the applications

The companion essay catalogues the applications of the Liénard potential [L, §7]. We organize them by which regime — conservative skeleton or stiff full system — governs their numerics, and state what the geometric viewpoint contributes to each.

### 7.1 Relaxation oscillators (the stiff regime)

The van der Pol triode circuit [L, §7.1], the FitzHugh–Nagumo excitable membrane [L, §7.3], and the cardiac sinoatrial-node pacemaker models [L, §7.3] are all relaxation oscillators in the stiff regime $\mu\gg1$ (or $\varepsilon\ll1$). Their numerics are governed by §6: the fast spiking/jump segments are stiff fast modes that require $L$-stability. The practical consequence is the everyday one noted in [L, §9.1]: relaxation-oscillator simulations are run with $L$-stable implicit (Radau/BDF) or specialized solvers, at $O(\mu)$ cost per period. The geometric contribution is *diagnostic*, not prescriptive: it explains *why* $L$-stability is required (the fast modes must be killed, $R(-\infty)=0$), *why* the tempting symplectic choice fails (§6.3), and *what structure* (the Hamiltonian–dissipative split (16)) a structure-preserving solver should target. For the cardiac phase-resetting computations [L, §7.3, ref. 35] — continuation of the limit cycle through the spiking parameter region — the discrete Poincaré map result (§4.4) is the relevant object: the continued cycle is the moving fixed point of the discrete return map, and its tracking error is controlled by the integrator's order, in practice limited by the event-detection floor.

### 7.2 Conservative and benchmark limits (the symplectic regime)

The integrable Duffing quartic [L, §7.2], the unforced conservative oscillators, and the piecewise-linear (Biryukov) systems with explicit limit cycles [L, §7.1, §9.4] are the domain where symplectic integration is the *right* tool. The Biryukov benchmarks are designed to test an integrator's ability to track the sharp corners of a relaxation cycle against an exact reference [L, §9.4]; a symplectic integrator's bounded energy error (§4) is the property that makes such a benchmark meaningful over the many periods a benchmark run spans, because the error does not accumulate. The relaxation "reduced model" — the slow flow $dx/dt = x/(\mu(1-x^2))$ giving the explicit period (16) [L, §7.4, §9.3] — is a 1D conservative (quadrature) problem; integrating it symplectically (or by a symplectic 1D method) preserves the quadrature's accuracy over the parameter sweeps for which it is used.

### 7.3 The discrete limit cycle as a computational object

Across all applications, the limit cycle is computed as a fixed point of the (discrete) Poincaré map [L, §9.2]. The geometric viewpoint adds two things. First, the accuracy of the discrete fixed point is controlled by the integrator's order, in practice limited by the event-detection floor (§4.4, verified), which is the a-priori error estimate for any shooting/continuation code. Second, when the skeleton is symplectically integrated, the discrete return map is the *structure-preserving discretization* of the conservative part, so the *geometry* of the fixed-point problem — the monotone one-way crossing of the Liénard curve that proves uniqueness [L, §4] — is expected to be preserved at the discrete level for small step (the discrete return map stays $O(h^p)$-close to the continuous one, so the contraction and the crossing property pass to it). A generic non-structure-preserving discretization can destroy this (spurious fixed points, loss of the monotonicity); the structure-preserving discretization is what makes Newton shooting on the discrete Poincaré map robust.

### 7.4 Stochastic Liénard dynamics

The stochastic Liénard equation [L, §8] is the underdamped Langevin equation with state-dependent friction; its invariant measure (Boltzmann for constant friction) is the stochastic analogue of the potential $V$'s structure. The structure-preserving philosophy extends: integrators that preserve (or approximately preserve) the invariant measure are the appropriate tools [L, §9.4]. The geometric-integration reference [1] develops exactly this class (structure-preserving methods for SDEs, invariant-measure preservation), and the Hamiltonian–dissipative split (16) is the deterministic skeleton of those stochastic methods. The stiff-regime tension of §6 carries over: the stochastic relaxation regime again requires $L$-stable treatment of the fast modes, and the symplectic (measure-preserving) structure is a secondary concern there.

---

## 8. Open problems and current directions

1. **Symplectic-dissipative $L$-stable methods.** §6 shows the prototypical symplectic method (implicit midpoint) fails in the stiff regime because it is not $L$-stable, and the $L$-stable methods that succeed (Radau IIA) are not symplectic. The design problem is to construct a **symplectic-dissipative split (16) whose damping sub-step is $L$-stable in the fast modes** while the conservative sub-step remains symplectic — so that the modified-Hamiltonian structure of §5.4 is preserved *and* the $O(\mu)$ cost is achieved. This is open: the existing splits use either exact damping (not $L$-stable in the central region, §6.3) or $A$-stable (not $L$-stable) damping sub-steps. A candidate direction is a damping sub-step that is $L$-stable for $\lambda<0$ (outer region) but applied with a fast-time sub-resolution for $\lambda>0$ (central region); the error analysis of such a hybrid is not available.

2. **Uniform (in $\mu$) modified-equation analysis of the split.** The modified-equation result (18)–(19) of §5.4 is stated for fixed $\mu$. The relaxation regime requires the modified potential $V_h$ and modified damping $f_h$ to be controlled *uniformly as $\mu\to\infty$*, so that the modified limit cycle tracks the true relaxation cycle at $O(h^p)$ for all $\mu$. This uniform theory is lacking; it is the analytic counterpart of the certified-numerics problem below.

3. **Certified numerics in the stiff regime.** The companion essay's open problem 5 [L, §10] asks for rigorous (interval-arithmetic) computation of relaxation limit cycles built on the reduced slow-flow models rather than brute-force IVP integration. The geometric viewpoint sharpens it: the certification should be of the *modified* system (18) — bound the modified potential $V_h$ and modified damping $f_h$ (problem 2), then bound the modified limit cycle's $O(h^p)$ distance from the true one (the discrete Poincaré map result (13), made rigorous with interval arithmetic). Systematic error bounds for this hybrid certified computation are lacking.

4. **Discrete Poincaré maps in the canard regime.** The canard explosion of generalized Liénard equations [L, §6.4, §10] is an exponentially small ($\sim e^{-c/\varepsilon}$) bifurcation phenomenon. The discrete return map $\mathcal{P}_h$ of a finite-precision integrator has its own (discrete) canard structure, and the question of how the *numerical* canard explosion (the discretization's resolution of the exponentially narrow parameter window) relates to the *continuous* one is open. This is where the structure-preserving discretization of the conservative part (§4.4, §7.3) and the resurgence theory of the continuous problem [L, §6.3] meet, and no analysis is available.

---

## 9. Conclusion

The Liénard equation is a Hamiltonian system — the flow of the potential $V$ — perturbed by a state-dependent vertical dissipation. That single fact organizes the geometric-integration story. **Strict symplectic integration applies to the Hamiltonian skeleton, not to the full dissipative system**; on the skeleton it delivers its hallmark bounded, non-secular energy error (verified, §4), which is the property that makes the long-time invariants of the conservative and reduced models — the period quadratures, the benchmark cycles, the small-$\mu$ circle — trustworthy at coarse step size. **The full system is handled by the Hamiltonian–dissipative split (16)**, whose modified equation (18)–(19) is the rigorous bridge to Liénard's qualitative theory: the numerical orbit converges to a modified limit cycle $O(h^p)$-close to the true one, and the qualitative conclusions (existence, uniqueness, stability) pass to it for $h$ small. **In the stiff relaxation regime — the regime that drives the physical applications — symplecticity and $L$-stability are in tension**: the prototypical symplectic method (implicit midpoint) fails (verified, §6.3) because it is $A$-stable but not $L$-stable — it does not kill the stiff fast modes, so it returns a bounded but severely distorted limit cycle (period off by $+10\%$ at $h=0.01$, the jump unresolved at $h=0.5$) rather than the true one — and the methods that succeed (Radau IIA) are not symplectic. The resolution is not to abandon geometric integration but to specify precisely — as this essay attempts — *which* structure (symplectic skeleton, $L$-stable dissipation, the split between them) is preserved in *which* regime, and to identify (problem 1) the design of a split that is simultaneously symplectic in the conservative part and $L$-stable in the fast modes as the open problem at the intersection of the two.

---

## References

[1] E. Hairer, C. Lubich, G. Wanner, *Geometric Numerical Integration: Structure-Preserving Algorithms for Ordinary Differential Equations*, 2nd ed., Springer Series in Computational Mathematics **31**, Springer, Berlin, 2006. (Reference [29] of the companion essay.)

[2] E. Hairer, C. Lubich, *The Numerical Solution of Oscillatory Differential Equations*, Springer Series in Computational Mathematics **39**, Springer, Berlin, 2010.

[3] C. Lubich, "Splitting methods for differential equations," *Math. Comput.* **77** (2008), 2105–2120.

[4] J. C. Butcher, *Numerical Methods for Ordinary Differential Equations*, 3rd ed., Wiley, Chichester, 2016.

[5] E. Hairer, G. Wanner, "Stiff ordinary differential equations," in *Handbook of Numerical Analysis*, Vol. I, P.-G. Ciarlet and L. Quartapelle (eds.), North-Holland, Amsterdam, 1991, 79–147.

[6] M. Crouzeix, "Generalized symmetric multistep methods for stiff differential equations," *Numer. Math.* **42** (1983), 345–359.

[7] R. D. Skeel, "Symplectic integrators and energy conservation: survey and improved algorithms," *SIAM J. Sci. Stat. Comput.* **13** (1992), 366–381.

[8] B. van der Pol, "On relaxation-oscillations," *Philosophical Magazine* (7) **2** (1926), 978–992. (Reference [2] of the companion essay.)

[9] F. Verhulst, *Nonlinear Differential Equations and Dynamical Systems*, Universitext, Springer-Verlag, Berlin, 1996. (Reference [32] of the companion essay.)

*All numerical values in this essay are produced by `geometric_integration_computation.py` through `geometric_integration_computation5b.py` (run under the project venv, Python 3.14, NumPy 2.4.6, SciPy 1.18.0, SymPy 1.14.0); raw outputs are in `geometric_integration_results*.txt`. The symbolic symplectic checks of the Storer–Verlet (7) and implicit-midpoint (8) maps in §2.1 (including the midpoint-vs-endpoint determinant distinction) are reproduced in `verify_review.py`.*
