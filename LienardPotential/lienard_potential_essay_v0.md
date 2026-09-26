# Liénard's Potential and Its Applications

*An essay on the Liénard system from the viewpoints of mathematical physics, functional and numerical analysis*

---

## 1. Introduction and terminology

Consider the second-order ordinary differential equation

$$
\ddot{x} + f(x)\,\dot{x} + g(x) = 0 , \qquad (1)
$$

with $f, g \in C^1(\mathbb{R})$, $f$ even and $g$ odd. Equation (1) is the **Liénard equation**, named after the French physicist Alfred-Marie Liénard, who introduced it in his 1928 study of maintained (self-sustained) electrical oscillations [1]. It is the prototypical model of a **nonlinear self-oscillator**: a damped mechanical or electrical system in which the damping is state-dependent in such a way that small oscillations are amplified while large ones are damped out, producing an isolated periodic orbit — a *limit cycle* — in the phase plane.

The two scalar functions entering (1) play structurally different roles, and it is worth fixing the vocabulary. Write the restoring force as the gradient of a one-variable potential,

$$
g(x) = V'(x), \qquad V(x) := \int_0^x g(s)\,ds , \qquad (2)
$$

and set

$$
F(x) := \int_0^x f(s)\,ds . \qquad (3)
$$

In this essay the term **"Liénard's potential"** refers to $V$ in (2), the potential whose gradient supplies the conservative (restoring) force of the system; the antiderivative $F$ of the damping coefficient is, in the classical literature, the **Liénard curve**. (We note at once that this usage is distinct from the *Liénard–Wiechert potentials* of relativistic electrodynamics, with which it has no connection beyond the author's name.)

Why is $V$ the central object of the theory? Three reasons. First, $V$ supplies the Hamiltonian skeleton of the dynamics: with $y = \dot{x}$, the energy

$$
E(x,y) = V(x) + \frac{y^2}{2}
$$

satisfies, along solutions of (1),

$$
\frac{dE}{dt} = V'(x)\,y + y\,\ddot{x} = g(x)\,y + y\bigl(-f(x)\,y - g(x)\bigr) = -f(x)\,y^2 . \qquad (4)
$$

Thus the Liénard equation is a Hamiltonian system with potential $V$ subjected to a *state-dependent dissipation*: energy is injected where $f(x) < 0$ and dissipated where $f(x) > 0$. The equilibria, the small-oscillation frequencies, and the multistability structure are all read off from $V$; $F$ then determines how the dissipation is distributed in phase space. Second, the qualitative theory of limit cycles — existence, uniqueness, and in particular the *number* of limit cycles — is governed by the interplay of the shapes of $V$ and $F$ (Sections 4–6). Third, and this is the point most relevant to numerical analysis, the asymptotic regimes of the equation — in particular the *relaxation regime*, where the period of the limit cycle is proportional to a large parameter and the orbit consists of slow segments glued by fast jumps — are controlled by the curvature of $V$ together with the folds of the Liénard curve $F$ (Sections 7.4 and 9).

A disambiguation is useful before proceeding. In much of the literature one fixes the harmonic potential $V(x) = x^2/2$ (so $g(x)=x$) and allows the Liénard curve $F$ to vary; (1) is then called the **classical Liénard equation** and written in the equivalent first-order form, obtained with the Liénard substitution $y = \dot{x} + F(x)$,

$$
\dot{x} = y - F(x), \qquad \dot{y} = -g(x) . \qquad (5)
$$

Both forms are used in this essay.

**Roadmap.** Section 2 recalls the historical origins. Section 3 develops the mathematical formulation and the standard examples. Section 4 states Liénard's theorem — the classical existence-and-uniqueness result — and its period computations. Section 5 describes the Lins Neto–de Melo–Pugh conjecture on the number of limit cycles and its resolution. Section 6 gives the functional-analytic perspective: the Poincaré map, Abelian integrals, the Chebyshev method, finiteness theorems, and the slow–fast (canard) formalism. Section 7 collects the applications in mathematical physics, including the relaxation limit with its explicitly computed period. Section 8 discusses the stochastic theory. Section 9 treats numerical analysis: stiffness, computation of limit cycles, and structure-preserving discretizations. Section 10 lists open problems.

---

## 2. Historical origins (1877–1928)

The equation descends from the acoustics of Lord Rayleigh. In his 1877 paper "On the maintenance of vibrations by forces of double period" [4] and, more systematically, in "On maintained vibrations" (1883) [5], Rayleigh modeled the sound of reed instruments by a linear oscillator with damping depending on the square of the velocity, and derived — by the method that now bears his name in the form of *Rayleigh's equation* $\ddot{x} + \delta(\dot{x}^2 - 1)\dot{x} + x = 0$ — the first rigorous analysis of an isolated periodic solution of a nonlinear ODE.

Georg Duffing, in his 1918 monograph on forced oscillations with amplitude-dependent natural frequency [3], studied the equation $\ddot{x} + \delta \dot{x} + \alpha x + \beta x^3 = \gamma \cos \omega t$: with $\alpha < 0 < \beta$ the potential $V(x) = \alpha x^2/2 + \beta x^4/4$ has two wells, and the system displays the coexistence of multiple stable states that makes the Duffing oscillator a model for buckled beams, magnetic pendula, and later for chaotic dynamics.

The direct line to the Liénard equation runs through radio physics. Liénard's 1928 paper "Étude des oscillations entretenues" [1] analyzed self-sustained oscillations in discharge circuits (and in the triode circuits being developed at the time), isolating precisely the structure (1) with even $f$ and odd $g$, proving the existence-and-uniqueness theorem of Section 4, and introducing the phase-plane transformation (5) with the Liénard curve $F$. Van der Pol's 1926 paper "On relaxation-oscillations" [2], derived from a mathematical model of a triode circuit, studied the case $\mu \gg 1$ of

$$
\ddot{x} - \mu(1 - x^2)\dot{x} + x = 0 \qquad (6)
$$

and coined the notion of a **relaxation oscillation**: an oscillation whose period grows proportionally to the nonlinearity parameter and whose waveform is a succession of slow drifts and rapid jumps. The vacuum-tube era (roughly 1920–1950) made (1) an engineering workhorse; its systematic dynamical-systems treatment — existence, cyclicity, bifurcations — was developed in parallel in the Western tradition (Coppel; Perko [8]) and in the Russian school (Andronov and Leontovich, and later Ye Yanqian and collaborators [23]).

Two programmatic facts framed the subsequent theory. Hilbert's sixteenth problem (1900) [34] asked for the maximal number $H(n)$ of limit cycles of a planar polynomial vector field of degree $n$, and their relative positions: the Liénard systems are a natural test class. And Dulac's problem, whether a *graphic* (a closed trajectory made of finitely many singularities and separatrices) can carry infinitely many limit cycles, was answered negatively by Ilyashenko [21] and Écalle [22] at the end of the 1980s, proving that limit cycles of analytic planar fields are finite in number.

---

## 3. Mathematical formulation

**Definition.** A *Liénard equation* is an equation of the form (1) with $f, g \in C^1(\mathbb{R})$, $f$ even, $g$ odd, and $x\,g(x) > 0$ for $x \neq 0$. With $V(x) = \int_0^x g(s)\,ds$ and $F(x) = \int_0^x f(s)\,ds$, the equivalent planar **Liénard system** is (5).

The equivalence is the Liénard transformation. Set $y = \dot{x} + F(x)$; then $\dot{y} = \ddot{x} + f(x)\,\dot{x} = -g(x)$, and $\dot{x} = y - F(x)$; conversely, given a solution of (5), $x(t)$ solves (1). The transformation is a near-identity change of variables in velocity: geometrically, the Liénard plane $(x, y)$ differs from the physical plane $(x, \dot{x})$ by the vertical shift $\dot{x} = y - F(x)$; the curve $y = F(x)$, the *Liénard curve*, is the $x$-nullcline of the Liénard system, and its folds are the turning points that organize the relaxation dynamics (Section 7.4).

Three structural observations are worth recording.

**(i) Energy balance.** Identity (4), $\dot E = -f(x)\dot x^2$, is the fundamental variational statement. The potential $V$ fixes a conserved quantity for the undamped skeleton ($f \equiv 0$), whose level sets $E = h$ are closed ovals around each nondegenerate minimum of $V$; the sign pattern of $f$ selects which annuli grow and which shrink. For (6) one has $f(x) = -\mu(1-x^2)$, so $\dot E = \mu(1-x^2)\dot x^2$: energy is injected on $|x| < 1$ and extracted on $|x| > 1$. No solution can stay in either region forever, and the limit cycle lives in between.

**(ii) Divergence.** The divergence of (5) is $\partial_x(y - F(x)) + \partial_y(-g(x)) = -F'(x) = -f(x)$, independent of $y$. Periodic orbits therefore feel the sign of $f$ directly: for a periodic orbit $\gamma$ of period $T$, the integral $\sigma = \int_\gamma \mathrm{div}\,dt$ determines the multiplier (in particular $\sigma < 0$ implies stability), which is the engine of most uniqueness proofs.

**(iii) Relation with Abel equations.** Eliminating $t$ with $v = v(x) = \dot{x}$ reduces (1) to

$$
v\,\frac{dv}{dx} + f(x)\,v + g(x) = 0 , \qquad (7)
$$

an Abel equation of the second kind. This is why explicit solutions are exceptional: first-order scalar equations of Abel type generally do not integrate in closed form, and the phase-plane (qualitative) methods of Sections 4–6 exist precisely because the direct method fails.

**Standard examples.**

- **Van der Pol oscillator** (6): $V(x) = x^2/2$, $F(x) = \mu(x^3/3 - x)$, $\mu > 0$. The Liénard curve is the cubic with folds at $x = \pm 1$. For every $\mu > 0$ there is a unique stable limit cycle; for $\mu \to 0$ it is a circle of radius 2 and period $2\pi + o(1)$, and for $\mu \to \infty$ it becomes a relaxation oscillation (Section 7.4).
- **Duffing oscillator** (unforced): $f \equiv \delta$, $V(x) = \alpha x^2/2 + \beta x^4/4$. With $\alpha < 0 < \beta$ the potential is double-welled and the origin is unstable, flanked by two stable equilibria; a single well ($\alpha > 0$) gives the usual hardening/softening oscillator. The forced equation exhibits subharmonic resonance and, for suitable parameters, chaotic motion [27].
- **FitzHugh–Nagumo system** (biological, Section 7.3):
  $$
  \dot{x} = \tfrac{1}{3}x^3 - x - y + b, \qquad \dot{y} = \varepsilon\,(x - a - y), \qquad 0 < \varepsilon \ll 1 . \qquad (8)
  $$
  With $Y = y - b$ and slow time $s = \varepsilon t$, (8) becomes
  $$
  \frac{dx}{ds} = -\frac{1}{\varepsilon}\bigl(Y - (x^3/3 - x)\bigr), \qquad \frac{dY}{ds} = -\bigl(x - (a+b) + Y\bigr) ,
  $$
  a slow–fast system of **(generalized) Liénard type**: the fast variable $x$ is slaved to the cubic nullcline $Y = x^3/3 - x$ — a cubic Liénard curve with folds at $x = \pm 1$ — while the slow variable drifts by $dY/ds = -(x - (a+b)) - Y$, the linear term $x - (a+b)$ playing the role of the restoring force and the $-Y$ term the (leaky) deviation from the classical Liénard form. In the singular limit $\varepsilon \to 0$ the trajectory is confined to the cubic, where the dynamics reduces to the scalar slow flow
  $$
  \frac{dx}{ds} = \frac{2x - (a+b) - x^3/3}{x^2 - 1} ,
  $$
  whose closed orbit (in the spiking parameter region) is the relaxation cycle of (8); the action potential is its fast segment.

---

## 4. Liénard's theorem: the potential as a criterion

The classical result, due to Liénard [1] and proved in modern form in [8, pp. 254–257], is:

**Theorem 1 (Liénard).** *Let $f, g \in C^1(\mathbb{R})$ with $f$ even, $g$ odd, $g(x) > 0$ for $x > 0$, and $F(x) = \int_0^x f(s)\,ds$ such that*

- *(i) $\lim_{x \to \infty} F(x) = +\infty$;*
- *(ii) $F$ has exactly one positive zero $p$;*
- *(iii) $F(x) < 0$ for $0 < x < p$, and $F(x) > 0$ and monotone increasing for $x > p$.*

*Then the Liénard system (5) has a unique limit cycle surrounding the origin, and it is stable (hyperbolic).*

The proof (Poincaré–Bendixson plus the Liénard-curve geometry) has two phases. For existence: the level curves of the energy $E = V + y^2/2$ rotate around the origin, and the sign pattern of condition (iii) makes the Liénard curve $y = F(x)$ a one-way membrane for the flow — a trajectory crossing it moves from the region where $f < 0$ (energy injection) into the region where $f > 0$ (energy extraction) — which confines the flow to an annulus; Poincaré–Bendixson then yields a periodic orbit. For uniqueness: one studies the return map of the flow to the branch of the Liénard curve with $x > p$; the geometric lemma — a trajectory leaving the curve at $x_0 > p$ re-enters it at $x_0' < x_0$, and the return map is strictly increasing, with a positive shift for small $x_0$ and a negative shift for large $x_0$ — forces its fixed point, the limit cycle, to be unique. The energy identity (4) is what makes the geometric lemma work, since the sign of $\dot E$ is prescribed pointwise by $f(x)$; the full argument is reproduced in [8, 9].

The role of the potential is visible in the hypotheses: $g(x) > 0$ for $x>0$ (with $g$ odd) is exactly the requirement that $V$ have a **single nondegenerate well** at the origin, $V(0) = 0$ the unique minimum; conditions (i)–(iii) then constrain $F$ to have a *single negative hump* followed by monotone escape to $+\infty$. In the van der Pol case, $F(x) = \mu(x^3/3 - x)$ has $p = \sqrt{3}$ and satisfies all three conditions, so (6) has a unique stable limit cycle for every $\mu > 0$.

A classical refinement is Sansone's period computation [10], with Massera's extension [11]: under the hypotheses of Theorem 1 the period $T$ of the unique limit cycle is bounded between two explicitly computable integrals in $F$ and $g$ (the period is the sum of quadratures over the two halves of the cycle, each half being the graph of a branch of the Abel equation (7) over a definite $x$-interval). This is one of the few cases where a *quantitative* statement about the limit cycle is available in closed quadrature — a fact of direct value to numerical analysis, since it provides an a priori check on computed periods.

**Remark.** The hypotheses of Theorem 1 are sufficient, not necessary; Liénard-type equations with asymmetric $f$ and $g$, or with several zeros of $F$, may still have a unique limit cycle, and the modern literature (see e.g. [9]) provides many weaker sufficient conditions. The theorem should be read as a *potential-shaped* criterion: single-well $V$ plus single-hump $F$ implies one stable cycle.

---

## 5. How many limit cycles? The Lins Neto–de Melo–Pugh conjecture

Fix the harmonic potential $V(x) = x^2/2$ and let the Liénard curve be a real polynomial of degree $n$:

$$
\dot{x} = y - F(x), \qquad \dot{y} = -x, \qquad F(x) = a_n x^n + \cdots + a_0 . \qquad (9)
$$

The question — how many limit cycles can (9) have? — was settled at the level of a conjecture by Lins Neto, de Melo and Pugh [12], who proved the sharp *lower* bound and conjectured the matching *upper* bound.

**Theorem 2 (Lins Neto–de Melo–Pugh, 1977 [12]).** *For every $n \ge 1$ there exist systems (9) of degree $n$ having at least $\lfloor (n-1)/2 \rfloor$ limit cycles.*

**Conjecture (Lins Neto–de Melo–Pugh).** *System (9) has at most $\lfloor (n-1)/2 \rfloor$ limit cycles.*

The lower bound is obtained by a clean averaging argument that is worth reproducing, since it is the archetype of the functional-analytic method of Section 6. Consider $\dot{x} = y + \varepsilon F(x)$, $\dot{y} = -x$ (a small perturbation of the center), in polar coordinates $x = r\cos\theta$, $y = r\sin\theta$:

$$
\frac{dr}{d\theta} = -\varepsilon \cos\theta\, F(r\cos\theta) + O(\varepsilon^2) .
$$

First-order averaging gives

$$
f(r) := \frac{1}{2\pi}\int_0^{2\pi} \cos\theta\, F(r\cos\theta)\,d\theta
= \sum_{j=0}^{\lfloor (n-1)/2 \rfloor} a_{2j+1}\,b_{2j+1}\, r^{2j+1} ,
$$

where only the odd part of $F$ survives and $b_{2j+1} = \frac{1}{2\pi}\int_0^{2\pi}\cos^{2j+2}\theta\,d\theta \neq 0$. Choosing the odd coefficients $a_1, a_3, \dots$ so that $f(r)$ has $\lfloor (n-1)/2 \rfloor$ simple positive zeros, the averaging theorem produces that many hyperbolic limit cycles, one near each zero. (This computation is exactly the leading-order term of the Poincaré–Pontryagin function of Section 6.2.)

The status of the conjecture, as surveyed in [13], is:

**Theorem 3 (status, per [13]).** *For system (9):*

- *(a) for $n = 1, 2$ there are no limit cycles;*
- *(b) for $n = 3, 4$ there is at most one limit cycle, and one exists;*
- *(c) for every $n \ge 6$ there exist systems (9) having at least $n - 2$ limit cycles.*

Parts (a) and (b) for $n=3$ were proved by Lins Neto–de Melo–Pugh [12] (divergence integrals plus the Poincaré–Bendixson annulus theorem); the case $n = 4$ — the uniqueness of the limit cycle for quartic $F$ — resisted for thirty-five years and was settled by Li and Llibre [17]. Part (c) is the story of the counterexamples: Dumortier, Panazzolo and Roussarie [14] first showed that the conjecture fails for $n \ge 7$ (producing $\lfloor(n-1)/2\rfloor + 1$ cycles); De Maesschalck and Dumortier [15] extended the failure to $n \ge 6$ with $\lfloor(n-1)/2\rfloor + 2$ cycles; and De Maesschalck and Huzak [16] finally proved the sharp-in-spirit lower bound $n-2$ for all $n \ge 6$, using *slow divergence integrals* (Section 6.4). The construction is by induction on the degree: one writes $F = F_{\mathrm{e}} + \delta F_{\mathrm{o}}$ (even plus a small odd part) and shows the slow divergence integral has an asymptotic expansion whose leading term is a polynomial with as many prescribed simple zeros as the induction requires; each simple zero then pins a limit cycle to the corresponding slow–fast cycle via the entry–exit relation, and the slow detuning $\lambda(\varepsilon)$ contributes one further cycle — $k$ zeros give $k+1$ cycles (Theorem 4). The even degrees $n \ge 6$ are built directly; the odd degrees $n \ge 7$ follow by perturbing the even-degree construction with a small additional odd term and applying the Poincaré–Bendixson annulus theorem to gain one more cycle.

Two consequences deserve emphasis. First, the conjecture as stated is **false**: for $n \ge 6$ the degree-$n$ Liénard systems admit essentially twice as many limit cycles as conjectured. Second, the case **$n = 5$ remains open** [13], as does the determination of the maximum number of limit cycles of (9) for every $n \ge 5$. Note the contrast with the role of the potential: all of this concerns systems with the *harmonic* potential. When $V$ is allowed to be a general even polynomial, the problem merges with Hilbert's sixteenth problem: Caubergh and Dumortier [20] studied precisely the classical Liénard equations of *even* degree from that viewpoint, and Ilyashenko and Panov [24] obtained explicit upper estimates for the number of limit cycles of Liénard equations in terms of the degrees of $F$ and $V$.

---

## 6. A functional-analytic perspective

### 6.1 The Poincaré map and the displacement function

Let $\Sigma$ be a transverse cross-section (say, the positive $x$-axis in the Liénard plane, or $x = 0$, $y > 0$). The **Poincaré (return) map** $\mathcal{P}: \Sigma \to \Sigma$ sends a point to its next intersection with $\Sigma$; limit cycles in a period annulus correspond bijectively to fixed points of $\mathcal{P}$, and their stability to the sign of $1 - \mathcal{P}'$. For the Liénard system, the *displacement function* $\Delta(z) = \mathcal{P}(z) - z$ is analytic on $\Sigma$ (for analytic $F, g$), and its zeros are isolated unless identically zero. The whole small-amplitude limit-cycle problem is thus a problem about the zeros of an analytic function defined by a flow — and the analytic structure of $\Delta$ is what the next two subsections exploit.

### 6.2 Abelian integrals and the Poincaré–Pontryagin function

Take the general potential $V$, analytic, with $V(0) = 0$, $V'(0) = 0$, $V''(0) > 0$, and consider the perturbation of the center

$$
\dot{x} = y - \varepsilon F(x), \qquad \dot{y} = -V'(x) , \qquad H(x,y) := \frac{y^2}{2} + V(x) . \qquad (10)
$$

The ovals $\mathcal{C}_h = \{H = h\}$, $0 < h < h_0$, fill a period annulus. The system (10) is the Hamiltonian system $(\dot x, \dot y) = (H_y, -H_x)$ perturbed by $(P, Q) = (-\varepsilon F(x), 0)$, and the **Poincaré–Pontryagin (Melnikov) function**

$$
I(h) := -\oint_{\mathcal{C}_h} F(x)\,dy , \qquad (11)
$$

is the first-order term of the displacement function: for every simple zero $h_i$ of $I$ there exists, for $\varepsilon$ small, a hyperbolic limit cycle of (10) near $\mathcal{C}_{h_i}$ (implicit function theorem applied to the exact displacement, which equals $\varepsilon I(h) + O(\varepsilon^2)$).

For the harmonic potential $V(x) = x^2/2$ the ovals are circles $x^2 + y^2 = 2h$, and (11) collapses to an elementary integral. Parameterizing $\mathcal{C}_h$ by $x = \sqrt{2h}\cos\theta$, $y = \sqrt{2h}\sin\theta$ and writing $F(x) = \sum a_j x^j$:

$$
I(h) = -\sqrt{2h}\sum_{k \ge 0} a_{2k+1}\,(2h)^k\, W_{2k+2},
\qquad
W_{2k+2} := \int_0^{2\pi} \cos^{2k+2}\theta\,d\theta
= \frac{2\pi\,(2k+2)!}{2^{2k+2}\,((k+1)!)^2} > 0 . \qquad (12)
$$

Only the odd part of $F$ contributes, and $I(h)$ is $\sqrt{2h}$ times a polynomial in $h$ of degree $\lfloor(n-1)/2\rfloor$. *Example.* For the van der Pol equation, $F(x) = -\mu(x - x^3/3)$, so $a_1 = -\mu$, $a_3 = \mu/3$ (an overall rescaling of $F$ only rescales $I$ and leaves its zeros fixed), and with $W_2 = \pi$, $W_4 = 3\pi/4$:

$$
I(h) = \mu\pi\sqrt{2h}\left(1 - \frac{h}{2}\right) ,
$$

whose unique positive zero is $h = 2$, i.e. a limit cycle on the circle $x^2 + y^2 = 4$ in the Liénard plane — radius 2. This is the classical result (limit cycle of amplitude 2 and period $2\pi + o(1)$ as $\mu \to 0$), and it is a good test of the bookkeeping: the same value $r = 2$ follows from first-order averaging (Section 5) and from direct numerical integration (Section 9.1).

For a general potential $V$, by contrast, (11) is a genuine **Abelian integral**: an integral over the oval $\mathcal{C}_h$ of a rational (or, for algebraic $V$, algebraic) one-form, depending analytically on the parameter $h$. Counting its isolated zeros is, in the Liénard case, the *Rokhlin problem* — the bounding of the number of zeros of the displacement function of an analytic perturbation of a center. The standard tool is the **Chebyshev method**: if the family of Abelian integrals $\{I_0(h), \dots, I_N(h)\}$ spanning the relevant space of perturbations forms an *extended complete Chebyshev (EC) system* on $(0, \varepsilon_0)$ — any nontrivial linear combination has at most $N$ zeros — then the number of small-amplitude limit cycles is bounded by $\dim$. EC properties are known for the period annulus of centers with even analytic $V$ under nondegeneracy assumptions on $V$; a modern, explicit Chebyshev criterion for Abelian integrals is due to Grau, Mañasas and Villadelprat [19]; the small-amplitude count for polynomial Liénard equations was carried out systematically by Blows and Lloyd [18]; and the classical monograph of the Russian school, Ye Yanqian and collaborators [23], is the standard reference for the Abelian-integral treatment of (1) and its generalizations.

### 6.3 Finiteness and Hilbert's sixteenth problem

The Chebyshev method bounds the number of limit cycles *near a center*. The global question — finiteness of the total number of limit cycles of a planar analytic (in particular polynomial) vector field — is the content of the **Écalle–Ilyashenko finiteness theorem**, which resolved Dulac's conjecture by proving that no graphic of an analytic planar vector field carries infinitely many limit cycles [21, 22]. The proof is a triumph of functional analysis applied to ODEs: Écalle's theory of *resurgent functions* (flattening, alien derivatives, resurgence) and Ilyashenko's analysis of the first-return map near the graphic — its asymptotics as an iterated product of *Ilyashenko functions*, flat functions of the form $e^{-\lambda/\varepsilon}\,\varepsilon^k\log^{\ell}\varepsilon$ — together reduce the cyclicity of a graphic to the finite rank of such a product in the scale of resurgent functions. Consequences: $H(n) < \infty$ for every $n$ (second part of Hilbert's sixteenth problem, finiteness part) [34], while the sharp values and even reasonable explicit upper bounds for $H(n)$ remain open; Ilyashenko and Panov [24] give explicit estimates in the Liénard case, and Caubergh and Dumortier [20] address the even-degree classical Liénard subclass.

### 6.4 Slow–fast functional analysis: slow divergence integrals and canards

Write the classical Liénard system in the singularly perturbed form

$$
\dot{x} = y - F(x), \qquad \dot{y} = -\varepsilon\, x , \qquad 0 < \varepsilon \ll 1 . \qquad (13)
$$

Assume $F(0) = F'(0) = 0$ and $F'(x)/x > 0$ for $x \neq 0$ (i.e. $F'$ has the same sign as $x$; for instance $F(x) = x^4/4$). The critical manifold $y = F(x)$ is attracting in the fast direction (since $\partial_x(y - F(x)) = -1$), and the *reduced* (slow) dynamics on it is $x' = -x/F'(x)$ (prime $= d/d\tau$, $\tau = \varepsilon t$). For each $x_0 > 0$ let $L(x_0) < 0$ be defined by $F(L(x_0)) = F(x_0)$; the piecewise-smooth closed curve made of the slow arc on $y = F(x)$ for $x \in [L(x_0), x_0]$ and the fast layer orbit closing it is a **slow–fast cycle** $\Gamma_{x_0}$. The **slow divergence integral** associated to $\Gamma_{x_0}$ ([13, (17)], following [16]) is

$$
I(x_0) := \int_{x_0}^{L(x_0)} \frac{F'(s)^2}{s}\,ds , \qquad (14)
$$

(the orientation of the limits is conventional — only the zero set of $I$ enters Theorem 4). It is, literally, the slow-time integral of the divergence along the slow arc: from $dx/d\tau = -x/F'(x)$ one gets $d\tau = -(F'(s)/s)\,ds$, hence $I(x_0) = \int_{\Gamma^s_{x_0}} \mathrm{div}\,X\,d\tau$, with $\mathrm{div}\,X = -F'(x)$ the divergence of the slow–fast vector field and the integral taken along the slow arc in the direction of the slow flow.

The key result, due to De Maesschalck and Huzak [16, Thm. 2], is:

**Theorem 4 (De Maesschalck–Huzak).** *Under $F(0) = F'(0) = 0$ and $F'(x)/x > 0$, if $I(x)$ has exactly $k$ simple positive zeros, then there is a smooth function $\lambda(\varepsilon)$ with $\lambda(0) = 0$ such that the perturbed system*
$$
\dot{x} = y - F(x), \qquad \dot{y} = \varepsilon\,(\lambda(\varepsilon) - x) \qquad (15)
$$
*has exactly $k + 1$ hyperbolic limit cycles, for $\varepsilon > 0$ small.*

The parameter $\lambda(\varepsilon)$ plays the role of a *slow-time detuning*: by shifting the slow drift one selects which slow–fast cycles are closed (the entry–exit relation is modified by $\lambda$), and each simple zero of $I$ pins one limit cycle to a prescribed slow–fast cycle. This is the machinery behind the counterexamples of Section 5: writing $F = F_{\mathrm{e}} + \delta F_{\mathrm{o}}$, one shows $I(x) = 2\delta I_1(x) + O(\delta^2)$ with an explicitly computable leading term $I_1$, and constructs $F$ of any even degree $n \ge 6$ for which $I_1$ has $n-3$ simple positive zeros, hence $n-2$ limit cycles [13, §5].

The same framework, pushed into the regime where the slow–fast cycles *touch* the folds of the critical manifold, produces **canard cycles**: orbits that follow the attracting branch, cross the fold, continue along the repelling branch for an $O(1)$ fast-time interval, and jump away. Dumortier and Roussarie's analysis of canard cycles and their bifurcations in generalized Liénard equations [25] and of relaxation bifurcations in the plane [26] shows that the number and stability of cycles can change in the exponentially small *canard explosion* regime, where the relevant parameter scales like $\exp(-c/\varepsilon)$ — a phenomenon that is invisible to any finite-order asymptotic expansion and is a genuine functional-analytic (resurgence-type) effect.

---

## 7. Applications in mathematical physics

### 7.1 Electrical circuits

The van der Pol equation (6) is derived from a triode amplifier with a nonlinear plate-current characteristic: the circuit equations reduce, in the appropriate dimensionless form, to a linear oscillator driven by a cubic nonlinearity in the displacement, i.e. (6) with $F(x) = \mu(x^3/3 - x)$ [2]. For $\mu \gg 1$ the oscillation is a **relaxation oscillation**: the period grows like $\mu$ and the waveform is nearly square, which is why such circuits were historically used as astable multivibrators and clock generators.

A fact of direct numerical value: if the damping $f$ is taken *piecewise constant* (so that $F$ is piecewise linear — the so-called Biryukov-type equations), the limit cycle is explicitly integrable, arc by arc [33]. These exactly-solvable Liénard systems serve as benchmark problems against which the accuracy of general numerical integrators can be measured (see also Section 9.4).

### 7.2 Mechanics and engineering

The unforced Duffing oscillator $\ddot{x} + \delta\dot{x} + \alpha x + \beta x^3 = 0$ is a Liénard equation with constant damping and potential $V(x) = \alpha x^2/2 + \beta x^4/4$ [3]. The double-well case ($\alpha < 0 < \beta$) models a buckled beam or a magnetic oscillator: two stable equilibria separated by an inverted barrier, with a (typically unstable) equilibrium at the origin; for $\delta = 0$ it is the integrable quartic oscillator, whose double-well phase portrait carries the separatrix (homoclinic) loops of the saddle at the origin. Under periodic forcing the Duffing equation displays subharmonic resonances, period-doubling routes to chaos, and strange attractors; it is the standard toy model of deterministic chaos in mechanical systems [27].

Liénard-type equations (with state-dependent friction $f$ and linear or weakly nonlinear $g$) also appear in seismological models of stick–slip fault motion and in models of vocal-fold oscillation in phonation; in both cases the mechanism is the same as in (4): a region of negative effective damping bounded by a region of positive damping.

### 7.3 Biology: excitable membranes and cardiac pacemakers

FitzHugh's 1961 reduction of the Hodgkin–Huxley model to a two-dimensional excitable system [6], and Nagumo, Arimoto and Yoshizawa's circuit implementation of pulse propagation along an "active" transmission line [7], led to the FitzHugh–Nagumo system (8). As shown in Section 3, (8) is a slow–fast system of (generalized) Liénard type: a cubic nullcline $Y = x^3/3 - x$ (folds at $x = \pm 1$) playing the role of the Liénard curve, and a slow drift $-(x - (a+b))$ playing the role of the restoring force. Its limit cycle, in the spiking parameter region, is a relaxation cycle whose fast segments model action potentials and whose slow segments model the recovery of the membrane; the singular limit $\varepsilon \to 0$ gives the piecewise-linear "geometric" action potential used in fast neural simulations. Cardiac pacemaker models in the same Liénard family (two-variable models of the sinoatrial node) are the setting of numerical phase-response computations such as [35], where the discontinuities in the phase-resetting response of a cardiac pacemaker are tracked by continuation methods — a concrete instance of Liénard dynamics meeting numerical bifurcation analysis.

### 7.4 The relaxation limit: an explicit singular-perturbation computation

The asymptotics of (6) as $\mu \to \infty$ is the cleanest explicit instance of the potential–dissipation interplay, and it is worth doing in full, since every step is checkable. Write (6) as $\dot{x} = y$, $\dot{y} = \mu(1 - x^2)y - x$. The variable $y$ is fast (its derivative carries $\mu$); the slow manifold is the zero set of the right-hand side in $y$,

$$
y_s(x) = \frac{x}{\mu(1 - x^2)} , \qquad |x| > 1 ,
$$

which is attracting because $\partial_y[\mu(1-x^2)y - x] = \mu(1 - x^2) < 0$ there. On it, the slow flow is

$$
\dot{x} = y_s(x) = \frac{x}{\mu(1 - x^2)} < 0 \quad (x > 1),
$$

so the right branch of the cycle is traversed from $x = 2$ down to the fold at $x = 1$ in time

$$
t_R = \mu \int_1^2 \frac{1 - x^2}{x}\,dx = \mu\left[\frac{x^2}{2} - \ln x\right]_1^2 = \mu\left(\frac{3}{2} - \ln 2\right) ,
$$

and by the symmetry $(x, y) \mapsto (-x, -y)$ the left branch takes the same time. The fast segments connect the lower end of the right slow branch ($x \approx 1$, where $y_s \to -\infty$) to the left branch at $x \approx -2$, and symmetrically $x \approx -1$ to $x \approx +2$: during the jump the trajectory crosses the central strip $|x| < 1$, where the fast direction is unstable ($\mu(1-x^2) > 0$), so the velocity $y$ is exponentially amplified before the trajectory is captured by the attracting branch at $x \approx \mp 2$; the jump lasts $O(1)$ in the fast time $\tau = \mu t$ and contributes only $O(1)$ to the period. Hence

$$
\boxed{\; T(\mu) = \mu\,(3 - 2\ln 2) + O(1) , \qquad \text{amplitude} \to 2 \; } \qquad (16)
$$

as $\mu \to \infty$ (classical result of van der Pol [2]; standard treatments in [27, 32]).

These asymptotics were verified in this essay by direct numerical integration (implicit Radau, tolerances $10^{-10}$, period measured from successive zero-crossings of $x$; the $\mu = 100$ period from successive jump times):

| $\mu$ | $T$ (measured) | $T/\mu$ | $T/(\mu(3-2\ln 2))$ | amplitude |
|---:|---:|---:|---:|---:|
| 0.05 | 6.2842 | 125.68 | 77.88 | 2.0000 |
| 0.1 | 6.2871 | 62.87 | 38.96 | 2.0001 |
| 1 | 6.6633 | 6.663 | 4.129 | 2.0086 |
| 10 | 19.078 | 1.908 | 1.182 | 2.0143 |
| 50 | 82.508 | 1.650 | 1.023 | 2.0030 |
| 100 | 162.84 | 1.628 | 1.009 | 2.0013 |

with $3 - 2\ln 2 = 1.61371\ldots$. Three facts are visible. As $\mu \to 0$, $T \to 2\pi$ (measured $T(0.05) = 6.2842$ versus $2\pi = 6.28319$): the cycle is the circle of radius 2 of Section 6.2. As $\mu \to \infty$, $T/\mu$ converges to $3 - 2\ln 2$ from above, the $O(1)$ correction being $+1.82$ at $\mu = 50$ and $+1.47$ at $\mu = 100$. The amplitude hovers near 2 in both limits (maximum $2.014$ at $\mu \approx 10$). The slow branches of the $\mu = 100$ orbit lie on $y_s(x)$ to within $2 \times 10^{-7}$, and the fast jumps run from $x \approx +1$ to $x \approx -2$ in $0.194$ time units (i.e. $\approx 19$ fast-time units), with a peak velocity $\max |\dot x| \approx 1.34\,\mu$ — the exponential amplification in the central strip quantified.

The canard regime (Section 6.4) is the fine structure of this picture in generalized Liénard equations: when a family of slow–fast cycles approaches the folds, an $O(\varepsilon)$ fast-time segment — a canard — is inserted, and cycles of mixed stability appear in an exponentially narrow parameter window [25, 26].

---

## 8. The stochastic perspective

Adding white noise to (1) gives the **stochastic Liénard equation**

$$
\ddot{x} + f(x)\,\dot{x} + V'(x) = \sigma\,\xi(t) , \qquad (17)
$$

the underdamped Langevin equation with state-dependent (nonlinear) friction. In the phase variables $(x, v)$, $v = \dot{x}$, the associated Fokker–Planck (Kramers) equation for the density $\rho(x,v,t)$ is

$$
\partial_t \rho = -v\,\partial_x\rho + f(x)\,\rho + \bigl(f(x)\,v + V'(x)\bigr)\,\partial_v \rho + \frac{\sigma^2}{2}\,\partial_v^2 \rho , \qquad (18)
$$

the forward equation of the SDE $\dot x = v$, $dv = -(f(x)v + V'(x))\,dt + \sigma\,dW_t$; it is the subject of Kramers' 1940 study of Brownian motion in a force field [28]. The potential $V$ reappears in the stationary solution: in the *linear*-friction case $f \equiv \gamma$ the invariant density is Boltzmann, $\rho_{\mathrm{eq}} \propto \exp\bigl(-(v^2/2 + V(x))/D\bigr)$ with $D = \sigma^2/(2\gamma)$, so $V$ directly prescribes the metastable structure. Two physical problems then become potential problems. *Kramers' escape*: for a double-well $V$ the mean first-passage rate over the barrier is $\sim \exp(-\Delta V/D)$ (to leading order), the classical result of [28]; with $f(x) \neq \gamma$ the stationary density is no longer a closed-form Boltzmann measure and the exponential is controlled by a friction-dependent action that reduces to $\Delta V/D$ for constant friction. *The overdamped limit*: formally removing the inertia (the $\gamma \to \infty$ reduction, on the slow time scale $t/\gamma$ with the effective noise strength $\sigma^2/(2\gamma)$ held fixed) (18) reduces to the one-dimensional Fokker–Planck equation on $x$ with drift $-V'$, and the Liénard limit cycle — a purely inertial phenomenon — disappears from the reduced dynamics (a one-dimensional gradient-type dynamics admits no isolated periodic orbits), a cautionary example of how the potential alone does not determine the oscillatory behavior: it is the *pair* $(V, F)$ that does.

Noise also acts on the limit cycle itself: a stable Liénard cycle in (17) becomes a noise-driven oscillation with a well-defined phase, and its phase-diffusion coefficient and stochastic bifurcations (noise-induced destruction or creation of the cycle) can be analyzed by the standard perturbation theory for perturbed planar flows. In the cardiac and neural models of Section 7.3, which are intrinsically stochastic in their physiological setting, the deterministic Liénard limit cycle is the skeleton around which the stochastic dynamics is organized; the numerical computation of phase-resetting curves for such models [35] is where this theory meets the numerical analysis of Section 9.

---

## 9. Numerical analysis

### 9.1 The regime map and stiffness

The numerical behavior of a Liénard equation is organized by the size of the nonlinearity parameter. For (6) there are three regimes.

**(i) Small $\mu$.** The cycle is $O(\mu)$-close to the circle of radius 2 in the Liénard plane (Section 6.2); the problem is a mildly nonlinear periodic IVP with no stiffness. Any standard explicit Runge–Kutta method with a step size modest relative to the period $2\pi$ suffices.

**(ii) Moderate $\mu$.** The cycle is rounded relaxation-like; explicit methods still work, but the step must begin to respect the fastest local time scale.

**(iii) The relaxation regime, $\mu \gg 1$.** This is where numerical analysis bites. The fast subsystem has eigenvalue $\lambda(x) = \mu(1 - x^2)$, so $|\lambda|_{\max} = 3\mu$ (at $x = \pm 2$), while the period of interest is $T \sim 1.614\,\mu$ (16). The **stiffness ratio** — period divided by fast time scale — is therefore of order $\mu^2$. Concretely: explicit RK4 is stable only for $\Delta t \cdot 3\mu \lesssim 2.785$, i.e. $\Delta t \lesssim 0.93/\mu$, and integrating one period requires $\gtrsim T/\Delta t \approx 1.7\,\mu^2$ steps — two orders of magnitude more than the $O(\mu)$ steps suggested by the period itself. L-stable implicit methods (Radau IIA, BDF) remove the stability constraint: a step $\Delta t = O(1)$ resolves the slow drift accurately and *passes through* the fast jumps (which last $O(1)$ fast-time units, i.e. $O(1/\mu)$ in physical time) without resolving them pointwise, reducing the cost per period to $O(\mu)$. This is not a theoretical curiosity; it is the everyday reason that relaxation-oscillator simulations (neural spiking, cardiac pacing, multivibrators) are run with implicit or specialized solvers.

The data of Section 7.4 were produced in this regime with an implicit Radau solver (tolerances $10^{-10}$), and they simultaneously validate the asymptotic theory: the $O(\mu^2)$-cost regime is exactly the one in which the $O(\mu)$-period of (16) can be measured at all.

### 9.2 Computing the limit cycle

The standard workflow has three ingredients [27, 30, 31].

**Shooting with Newton on the Poincaré map.** Choose a cross-section $\Sigma$; the map $\mathcal{P}$ is computed by integrating one full orbit; a fixed point is found by Newton's method on $\mathcal{P}(z) - z$, with the $2 \times 2$ (or higher, in extended problems) Jacobian obtained by the variational equation. Each Newton iteration costs one orbit integration, and in the relaxation regime that orbit is stiff — the $O(\mu^2)$ cost of Section 9.1 is inherited by every iteration.

**Boundary-value (periodic normalization) formulation.** The orbit is posed as a periodic boundary-value problem on one period, with a normalization condition (e.g. fixing the phase by a linear constraint on the period) that makes the linearized operator invertible. Kuznetsov, Govaerts, Doedel and Dhooge [30] developed the periodic-normalization framework for codimension-one bifurcations of limit cycles; Beyn et al. [31] treat the general continuation and normal-form computations. The BVP formulation is advantageous precisely in the stiff regime, because the stiffness is distributed over the whole period and handled by a collocation or BDF-type discretization rather than by a single forward IVP.

**Continuation in parameters.** Limit cycles are continued in the nonlinearity (or a physical) parameter by pseudo-arc-length continuation, with local bifurcations (Hopf, homoclinic, period-doubling, canard points) detected from the eigenvalues of the linearized Poincaré map. This is the content of the AUTO family of codes and of the modern packages built on [30, 31].

### 9.3 Asymptotics as a computational tool

The two explicit asymptotics of Sections 5 and 7.4 double as numerical methods. The small-$\mu$ averaged equation (radius 2, period $2\pi$) provides the initial guess for Newton shooting at moderate $\mu$; the large-$\mu$ slow flow $dx/dt = x/(\mu(1-x^2))$ (16) gives, to $O(1/\mu)$ relative accuracy in the period, a closed-form trajectory that can be used as a *reduced model* for fast parameter sweeps. Verhulst's systematic treatment of relaxation oscillations, canards, and their asymptotic expansions [32] is the standard reference for the error estimates of such reductions (in particular the $O(1)$ period correction of (16) and its higher-order terms).

### 9.4 Structure-preserving discretization

A Liénard system is not Hamiltonian — $\mathrm{div} = -f \not\equiv 0$ — so symplectic integration, in the strict sense, does not apply; but the energy identity (4) still prescribes a natural structure to preserve: the splitting, in the physical variables $(x, v)$ with $v = \dot x$, between the conservative part $(\dot x, \dot v) = (v, -V'(x))$ and the dissipative part $(0, -f(x)\,v)$. Splitting-type schemes (integrate the Hamiltonian part with a symplectic step, treat the damping exactly or with an A-stable sub-step) inherit the long-time energy behavior of the exact flow and avoid the spurious amplitude drift that plagues generic high-order explicit methods on long runs. The general theory of such structure-preserving algorithms — symplectic, energy-stable, and invariant-measure-preserving methods, analyzed via modified equations and backward error analysis — is developed in [29]. For the stochastic equation (17) the same philosophy applies at the level of the invariant measure of the Kramers equation (18): integrators that preserve (or approximately preserve) the Boltzmann structure of the stationary density are the appropriate tools [29, 28].

Finally, the piecewise-linear (Biryukov-type) Liénard systems, whose limit cycles are explicit [33], are natural *benchmark problems*: they test an integrator's ability to track the sharp corners of a relaxation cycle, where the error of a method is maximally exposed, against a reference that is exact.

---

## 10. Open problems

1. **The case $n = 5$.** The Lins Neto–de Melo–Pugh conjecture holds for $n \le 4$ and fails for $n \ge 6$; whether a degree-five $F$ admits three limit cycles is open [13], as is the exact maximum for every $n \ge 5$.
2. **Hilbert's sixteenth problem.** The finiteness $H(n) < \infty$ is settled [21, 22], but the values (and relative positions) of the maximal numbers of limit cycles remain unknown for all $n \ge 2$ in the polynomial class; the Liénard subclass provides the sharpest estimates [20, 24].
3. **Cyclicity of graphics and the resurgence mechanism.** The Écalle–Ilyashenko theory proves finiteness but yields poor explicit bounds; understanding the *order* of cyclicity of nondegenerate graphics in concrete Liénard families — and the role of alien derivations in bounding it — is an active direction [21, 22].
4. **Canard explosions.** The exponentially small scaling of the bifurcation parameter in the canard explosion of generalized Liénard equations [25, 26] admits rigorous asymptotic descriptions in special cases; a general, uniform (in $\varepsilon$) theory of the explosion — including its detection and resolution in finite precision — is open.
5. **Numerical certification in the stiff regime.** The $O(\mu^2)$ cost of explicit methods (Section 9.1) suggests that rigorous numerics (validated computation with interval arithmetic) of relaxation limit cycles should be built on the reduced slow-flow models rather than on brute-force IVP integration; systematic error bounds for such hybrid certified computations are lacking.
6. **Stochastic Liénard dynamics.** Rigorous results on the phase diffusion of the noisy limit cycle (17) in the relaxation regime, and on invariant measures of the Kramers equation (18) with state-dependent friction, are largely perturbative; quantitative non-perturbative bounds are open.

---

## 11. Conclusion

The potential $V$ of a Liénard system is the object that carries the conservative content of the dynamics — its wells set the equilibria, its curvature sets the small-oscillation frequencies, and its ovals $\{H = h\}$ supply the annuli in which limit cycles live. The Liénard curve $F$ then decides where energy is injected and where it is withdrawn, and it is the interplay of the two that the entire theory organizes: Liénard's theorem reads single-well $V$ plus single-hump $F$ as a certificate of exactly one stable cycle [1, 8]; the Lins Neto–de Melo–Pugh conjecture and its partial resolution read the *degree* of $F$ (with $V$ harmonic) as a candidate for the maximal number of cycles [12–17]; the Abelian-integral and Chebyshev methods read the *period annulus of $V$* as the domain on which the displacement function lives [18, 19, 23]; the slow–fast theory reads the *folds of $F$ against the curvature of $V$* as the seats of canards and relaxation bifurcations [16, 25, 26]; and the numerical analysis reads the same geometry as a stiffness problem whose cost is $O(\mu^2)$ for explicit methods and $O(\mu)$ for implicit ones, and whose asymptotic limits are explicit quadratures [2, 27, 29, 32]. One equation, one potential, and — from 1877 to the present — a continuously deepening understanding of what the pair $(V, F)$ can do.

---

## References

[1] A. Liénard, "Étude des oscillations entretenues", *Revue Générale de l'Électricité* **23** (1928), 901–912, 946–954.

[2] B. van der Pol, "On relaxation-oscillations", *Philosophical Magazine* (7) **2** (1926), 978–992.

[3] G. Duffing, *Erzwungene Schwingungen bei veränderlicher Eigenfrequenz und ihre technische Bedeutung*, Vieweg, Braunschweig, 1918.

[4] J. W. Strutt (Lord Rayleigh), "On the maintenance of vibrations by forces of double period", *Philosophical Magazine* (5) **4** (1877).

[5] J. W. Strutt (Lord Rayleigh), "On maintained vibrations", *Philosophical Magazine* (5) **15**, no. 94 (1883). doi:10.1080/14786448308627342.

[6] R. FitzHugh, "Impulses and physiological states in theoretical models of nerve membrane", *Biophysical Journal* **1**(6) (1961), 445–466.

[7] J. Nagumo, S. Arimoto, S. Yoshizawa, "An active pulse transmission line simulating nerve axon", *Proceedings of the IRE* **50**(10) (1962), 2061–2070.

[8] L. Perko, *Differential Equations and Dynamical Systems*, 2nd ed., Texts in Applied Mathematics 7, Springer-Verlag, New York, 1991 (proof of Liénard's theorem: pp. 254–257).

[9] F. Dumortier, J. Llibre, J. C. Artés, *Qualitative Theory of Planar Differential Systems*, Universitext, Springer-Verlag, New York, 2006.

[10] G. Sansone, "Soluzioni periodiche dell'equazione di Liénard. Calcolo del periodo", *Rendiconti del Seminario Matematico dell'Università e del Politecnico di Torino* **10** (1951), 155–171.

[11] J. L. Massera, "Sur un théorème de G. Sansone sur l'équation de Liénard", *Bollettino dell'Unione Matematica Italiana* (3) **9** (1954), 367–369.

[12] A. Lins Neto, W. de Melo, C. C. Pugh, "On Liénard's equation", in J. Palis, M. do Carmo (eds.), *Geometry and Topology* (São Paulo, 1976), Lecture Notes in Mathematics **597**, Springer-Verlag, 1977, 335–357.

[13] J. Llibre, X. Zhang, "Limit cycles of the classical Liénard differential systems: A survey on the Lins Neto, de Melo and Pugh's conjecture", *Expositiones Mathematicae* **35**(3) (2017), 286–299. doi:10.1016/j.exmath.2016.12.001.

[14] F. Dumortier, D. Panazzolo, R. Roussarie, "More limit cycles than expected in Liénard equations", *Proceedings of the American Mathematical Society* **135** (2007), 1895–1904.

[15] P. De Maesschalck, F. Dumortier, "Classical Liénard equation of degree $n \ge 6$ can have $\lfloor (n-1)/2 \rfloor + 2$ limit cycles", *Journal of Differential Equations* **250** (2011), 2162–2176.

[16] P. De Maesschalck, R. Huzak, "Slow divergence integrals in classical Liénard equations near centers", *Journal of Dynamics and Differential Equations* **27**(1) (2015), 177–185.

[17] C. Li, J. Llibre, "Uniqueness of limit cycle for Liénard equations of degree four", *Journal of Differential Equations* **252** (2012), 3142–3162.

[18] T. R. Blows, N. G. Lloyd, "The number of small-amplitude limit cycles of Liénard equations", *Mathematical Proceedings of the Cambridge Philosophical Society* **95** (1984), 359–366.

[19] M. Grau, F. Mañasas, J. Villadelprat, "A Chebyshev criterion for Abelian integrals", *Transactions of the American Mathematical Society* **363** (2011), 109–129.

[20] M. Caubergh, F. Dumortier, "Hilbert's 16th problem for classical Liénard equations of even degree", *Journal of Differential Equations* **244** (2008), 1359–1394.

[21] Yu. S. Ilyashenko, *Finiteness Theorems for Limit Cycles*, Translations of Mathematical Monographs **94**, American Mathematical Society, Providence, RI, 1991.

[22] J. Écalle, *Introduction aux fonctions analysables et preuve constructive de la conjecture de Dulac*, Hermann, Paris, 1992.

[23] Ye Yanqian et al., *Theory of Limit Cycles*, Translations of Mathematical Monographs **66**, American Mathematical Society, Providence, RI, 1986.

[24] Yu. Ilyashenko, A. Panov, "Some upper estimates of the number of limit cycles of planar vector fields with applications to Liénard equations", *Moscow Mathematical Journal* **1** (2001), 583–599.

[25] F. Dumortier, R. Roussarie, "Multiple canard cycles in generalized Liénard equations", *Journal of Differential Equations* **174** (2001), 1–29.

[26] F. Dumortier, R. Roussarie, "Bifurcation of relaxation oscillations in dimension two", *Discrete and Continuous Dynamical Systems* **19** (2007), 631–674.

[27] Y. A. Kuznetsov, *Elements of Applied Bifurcation Theory*, 3rd ed., Applied Mathematical Sciences **112**, Springer, New York, 2004.

[28] H. A. Kramers, "Brownian motion in a field of force and the diffusion model of chemical reactions", *Physica (Utrecht)* **7** (1940), 284–360.

[29] E. Hairer, C. Lubich, G. Wanner, *Geometric Numerical Integration: Structure-Preserving Algorithms for Ordinary Differential Equations*, 2nd ed., Springer Series in Computational Mathematics **31**, Springer, Berlin, 2006.

[30] Y. A. Kuznetsov, W. Govaerts, E. J. Doedel, A. Dhooge, "Numerical periodic normalization for codim 1 bifurcations of limit cycles", *SIAM Journal on Numerical Analysis* **43**(4) (2005), 1407–1435.

[31] W. J. Beyn, A. Champneys, E. J. Doedel, W. Govaerts, Y. A. Kuznetsov, B. Sandstede, "Numerical continuation and computation of normal forms", in B. Fiedler (ed.), *Handbook of Dynamical Systems*, Vol. 2, Elsevier, 2002, 149–219.

[32] F. Verhulst, *Nonlinear Differential Equations and Dynamical Systems*, Universitext, Springer-Verlag, Berlin, 1996.

[33] A. M. Pilipenko, V. N. Biryukov, "Investigation of modern numerical analysis methods of self-oscillatory circuits efficiency", *Journal of Radio Electronics* **9** (2013).

[34] D. Hilbert, "Mathematische Probleme", *Göttinger Nachrichten* (1900), 253–297; English translation: *Bulletin of the American Mathematical Society* **8** (1902), 213–254.

[35] T. Krogh-Madsen, E. J. Doedel, L. Glass, M. R. Guevara, "Apparent discontinuities in the phase-resetting response of cardiac pacemakers", *Journal of Theoretical Biology* **230**(4) (2004), 499–519.
