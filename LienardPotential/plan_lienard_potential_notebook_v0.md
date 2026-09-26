# Plan — JupyterLab Notebook: *Liénard's Potential and Its Applications*

**Goal.** Generate, from `lienard_potential_essay_v3.md`, an executable JupyterLab
notebook that covers **every applicability area of the essay that is numerically or
symbolically treatable in Python**, each illustrated with an adequately
explanative, self-contained simulation whose output is *verified* (not assumed).
The notebook is a **computational companion** to the essay: it does not re-derive
the proof-theoretic material (canards, resurgence, the Hilbert-16th values), but it
makes every *computable* claim in the essay reproducible and visually legible.

**Deliverable.** A single `lienard_potential_notebook_v0.ipynb` that runs
end-to-end in the workspace `.venv` and opens in JupyterLab, plus the builder script
that deterministically produces it (the source of truth, per your part-file
convention).

---

## 0. Environment (verified)

| Item | Value | Status |
|---|---|---|
| Interpreter | `./.venv/Scripts/python.exe` (workspace root) | ✅ Python 3.14.3 |
| numpy | 2.4.6 | ✅ |
| scipy | 1.18.0 (`solve_ivp`: Radau, BDF, LSODA, RK45) | ✅ |
| sympy | 1.14.0 (symbolic integration, `solve`, `simplify`) | ✅ |
| matplotlib | 3.11.1 (Agg for headless exec) | ✅ |
| pandas | 3.0.5 (result tables) | ✅ |
| nbformat / nbconvert | 5.11.1 / 7.17.1 | ✅ |
| jupyterlab / jupyter | 4.6.3 / 1.1.1 | ✅ |
| mpmath | 1.3.0 (high-precision cross-checks) | ✅ |

**Solver strategy (fixed up front, used consistently):**

- **Non-stiff / moderate-μ regimes** → `RK45` or `DOP853` (explicit, cheap).
- **Relaxation / stiff regime (μ ≫ 1, FitzHugh–Nagumo, Biryukov)** → `Radau`
  (L-stable implicit). This is the whole point of §9.1 and is the *reason* the
  notebook uses Radau for the flagship computations.
- **Tolerances** → `rtol=atol=1e-10` for the headline van der Pol table (matches the
  essay's stated `10⁻¹⁰`); `1e-9`/`1e-11` elsewhere.
- **Period measurement** → median of successive *same-direction* (positive-going)
  zero-crossings of `x`, computed on the settled portion of the trajectory (discard
  transient). *Pitfall locked in from probing:* measuring between *any* crossings
  (or between same-direction crossings without discarding the first) returns a
  half-period — the first probe did exactly this and returned T/μ≈0.81 (half the
  correct 1.61). The notebook helper `period_from_crossings(t, x, direction=+1)`
  enforces the correct convention.
- **Runtime budget** → target **≤ 3 min total**. Heaviest cells: stochastic
  stationary sample (~1.5M Euler–Maruyama steps, ~10 s) and the 4-run Kramers sweep
  (~30–60 s). A `QUICK = False` flag at the top scales these down (fewer steps,
  coarser `dt`) so a first pass runs in ~30 s.

---

## 1. Scope: what is treatable, and how

The notebook covers the essay's applicability areas that reduce to (a) **integrating
an ODE/SDE** and/or (b) **a closed-form or quadrature/sympy computation**. The
proof-theoretic areas are *named and linked* but not simulated.

| Essay § | Topic | Treatable? | Method | Notebook § |
|---|---|---|---|---|
| 1, 3 | Liénard transform (1)↔(5), energy identity (4) | ✅ | **Symbolic** (sympy) | NB §2 |
| 3 | Standard examples (van der Pol, Duffing, FitzHugh–Nagumo) | ✅ | **Numeric** (phase portraits) | NB §3 |
| 4 | Liénard's theorem (single-well V + single-hump F ⇒ 1 stable cycle) | ✅ | **Numeric** (return map on a section) | NB §4 |
| 5 | Lins Neto–de Melo–Pugh: averaging lower bound | ✅ | **Symbolic** (averaging / Poincaré–Pontryagin) | NB §5 |
| 6.1–6.2 | Poincaré map, Abelian integrals, I(h) | ✅ | **Symbolic + numeric** (plot I(h), count zeros) | NB §5 |
| 6.4 | Slow divergence integral, Theorem 4 (machinery) | ✅ | **Numeric** (brentq + quad), demo only | NB §6 |
| 7.1 | Electrical circuits; Biryukov piecewise-linear benchmark | ✅ | **Numeric** (explicit cycle as benchmark) | NB §7 |
| 7.2 | Duffing double-well, separatrix, multistability | ✅ | **Numeric** (phase portrait, energy) | NB §8 |
| 7.3 | FitzHugh–Nagumo relaxation / action potentials | ✅ | **Numeric** (stiff, spiking params) | NB §9 |
| 7.4 | Relaxation limit: explicit period (16) | ✅ | **Numeric** (flagship Radau table) | NB §10 |
| 8 | Stochastic Liénard: Boltzmann density, Kramers escape | ✅ | **Numeric SDE** (Euler–Maruyama) | NB §11 |
| 9.1 | Stiffness / regime map: O(μ²) vs O(μ) | ✅ | **Numeric** (RK4 diverges vs Radau) | NB §10 |
| 9.4 | Structure-preserving (symplectic split) | ✅ | **Numeric** (leapfrog vs RK4 energy drift) | NB §12 |
| 2, 6.3 | History; Hilbert-16th values; Écalle–Ilyashenko finiteness | ❌ | Theoretical only — **prose cell, no sim** | NB §1, §6 |
| 6.4 (canard explosion) | Resurgence, exp(−c/ε) scaling | ❌ | Theoretical only — **prose cell** | NB §6 |
| 10 | Open problems (n=5, H(n), certified numerics) | ❌ | **Prose cell** (framing) | NB §13 |

> **Design rule.** Every numeric cell ends with a *verification line* that compares
> the computed quantity to a known analytic/asymptotic value and prints the residual.
> No cell is "trust me" — each either matches a closed form, an asymptotic limit,
> or a benchmark, and says so.

---

## 2. Notebook structure (cell-by-cell outline)

The notebook is assembled by a **builder script** (`build_lienard_notebook.py`,
written with `nbformat.v4`) so the content is versionable Python, not hand-edited
JSON. Sections below give, for each: the essay anchor, the code approach, the
**expected output (measured in the .venv during planning)**, and the figure.

### NB §0 — Setup & conventions
- Markdown: title, one-paragraph scope statement, pointer to the essay, the
  `(V, F)` vocabulary (V = Liénard's potential = gradient of restoring force;
  F = Liénard curve = antiderivative of damping).
- Code: `import numpy/scipy/sympy/matplotlib/pandas`; `matplotlib.rcParams`
  (Agg, 120 dpi, grid); the shared helpers:
  - `period_from_crossings(t, x, direction=+1)` — settled-window, same-direction.
  - `zero_crossings(t, x, level=0.0, direction)`.
  - `QUICK` flag controlling stochastic run lengths.
- **Expected:** imports succeed; helper self-test on a known sine (period = 2π/ω).

### NB §1 — The Liénard system in one page (§1, §3)
- Markdown: eqs (1)–(5); the two roles of V and F; energy identity (4)
  `dE/dt = −f(x)·ẋ²`; divergence `= −f(x)`.
- Code (symbolic, **NB §2 content folded in**): define generic `f, g, F`; use sympy
  to (a) verify the Liénard transform — `d/dt[x' + F(x)] = x'' + f(x)x' = −g(x)`
  after substituting `x'' = −f(x)x' − g(x)`; (b) verify `dE/dt = −f(x)ẋ²` for
  `E = V + ẋ²/2`.
- **Expected:** both sympy simplifications return the claimed identity (printed
  `True`). *Verified during planning* — the transform identity reduces exactly to
  `−g(x(t))`.

### NB §2 — Energy balance & the van der Pol skeleton (§3.i)
- Code: `f(x) = −μ(1−x²)`; plot the energy-injection/dissipation regions
  (shade `|x|<1` where `f<0`, `|x|>1` where `f>0`); overlay a short RK45 orbit and
  its `E(t)` to show `E` rises in the strip, falls outside.
- Figure: phase plane with shaded annulus + a decaying-inward / growing-outward
  orbit pair.
- **Expected:** `E(t)` clearly non-monotone, with net drift toward the cycle.

### NB §3 — Standard examples: three phase portraits (§3)
- Code: one figure, three panels:
  1. **van der Pol** μ=1 (unique stable cycle, radius ≈2).
  2. **Duffing** double-well (α=−1, β=1, δ=0) — saddle at 0, two wells, separatrix.
  3. **FitzHugh–Nagumo** (a=0, b=0, ε=0.01) — spiking relaxation cycle.
- Each panel: trajectory + equilibrium(s) marked + `F(x)` nullcline where relevant.
- **Expected (measured):** van der Pol max|x|≈2.01; Duffing equilibria 0,±1;
  FH-N spike height≈4, x∈[−1.99, 1.99].

### NB §4 — Liénard's theorem as a *certificate* (§4)
- Markdown: state Theorem 1 (single-well V + single-hump F ⇒ unique stable cycle);
  note it is sufficient, not necessary.
- Code: for van der Pol (F has one positive zero p=√3), build the **return (Poincaré)
  map** on the section `x=0, y>0`: integrate from a grid of starting points on the
  section until the next return; plot `P(z) − z` (the displacement) and mark its
  single zero = the limit cycle; confirm `P'(z*) < 1` (stability).
- **Expected:** exactly one sign change of the displacement → one cycle;
  multiplier < 1. *This is the numeric face of the theorem.*

### NB §5 — How many cycles? Averaging & the Poincaré–Pontryagin function (§5, §6.2)
- Code (symbolic, **the symbolic centerpiece**):
  - Compute `W_{2k+2} = ∫₀^{2π} cos^{2k+2}θ dθ` in sympy for k=0,1,2 → π, 3π/4,
    5π/8; cross-check against the closed form
    `2π(2k+2)!/(2^{2k+2}((k+1)!)²)`.
  - For van der Pol F (a₁=−μ, a₃=μ/3), build
    `I(h) = −√(2h) Σ a_{2k+1}(2h)^k W_{2k+2}` → simplify to
    `μπ√(2h)(1 − h/2)`; `solve(I(h)=0, h)` → `h=2` → radius `r=√(2h)=2`.
  - Generalize: for a *cubic* F with free odd coefficients, show `I(h)` is
    `√(2h)·(a₁ W₂ + a₃(2h) W₄)`, so its positive zeros (≤1 for cubic) = candidate
    small-amplitude cycles. Plot `I(h)` for a few coefficient choices and mark
    zeros.
- **Expected (measured):** W values and I(h) zero at h=2 / r=2 reproduce exactly.
  This is the first-order term of the displacement function and the archetype of
  the functional-analytic method.
- Markdown: connect to Theorem 2 (lower bound ⌊(n−1)/2⌋) and the (false)
  Lins Neto–de Melo–Pugh conjecture; note n=5 open.

### NB §6 — Slow–fast theory: the machinery, not the counterexample (§6.4)
- Code: implement the **slow divergence integral** `I(x₀) = ∫_{x₀}^{L(x₀)} F'(s)²/s ds`
  with `L(x₀)` found by `brentq` on `F(L)=F(x₀)` (L<0) and the integral by `quad`.
  Demo on a valid even F (`F = x⁴/4 + δx³`, δ=0.5): plot `I(x₀)`, mark its simple
  zeros, state "k zeros ⇒ Theorem 4 gives k+1 cycles" (here 1 zero ⇒ 2 cycles).
- Markdown: explain the entry–exit relation and the slow detuning λ(ε); **explicitly
  scope out** the full n−2-cycle counterexample construction and the canard-explosion
  (resurgence) regime — name them, cite [16,25,26], do not simulate.
- **Expected (measured):** integral evaluates; 1 zero found for the demo F.

### NB §7 — Biryukov piecewise-linear benchmark (§7.1, §9.4)
- Markdown: piecewise-constant damping ⇒ F piecewise linear ⇒ limit cycle
  **explicitly integrable arc by arc**; these are the benchmark problems of §9.4.
- Code: `f(x) = −1 (|x|<1)`, `+1 (|x|>1)`; integrate with Radau from *three
  different initial conditions*; show all three converge to the **same** cycle;
  report amplitude and period; overlay the piecewise-linear F.
- **Expected (measured):** IC-independent cycle, amplitude ≈ 2.5012, period ≈
  6.5658. *Verified to be IC-independent during planning* — this is the "exact
  reference" an integrator is measured against.

### NB §8 — Duffing double-well & multistability (§7.2)
- Code: `V = αx²/2 + βx⁴/4` (α=−1, β=1); plot V with wells at ±1 (V=−0.25) and the
  barrier at 0 (V=0); draw the **separatrix** (E=E_sep=0) homoclinic loops; show a
  small orbit in each well and the saddle at the origin; report equilibria and the
  separatrix energy.
- **Expected (measured):** equilibria 0,±1; V(±1)=−0.25; separatrix E=0;
  conservative orbit E-drift ~3e−10 (δ=0).
- Markdown: link to buckled beams / magnetic pendula; note forced Duffing → chaos
  (out of scope, prose only).

### NB §9 — FitzHugh–Nagumo: action potentials (§7.3)
- Markdown: eq (8); the cubic nullcline `Y = x − x³/3` (folds at ±1) as the Liénard
  curve; the spiking condition `x* = (3(a+b))^{1/3}` must lie on the middle branch
  (|x*|<1).
- Code: **parameter scan** — a grid of (a,b) colored by whether the orbit spikes
  (limit cycle) or settles (equilibrium); then a full time-series + phase portrait
  of a spiking case (a=0, b=0) with the action-potential spike highlighted.
- **Expected (measured):** spiking at (0,0) → period≈221, spike≈4; **non-spiking at
  (0.7, 0.8)** (x*=1.651 settles — the essay's "spiking region" is *not* all
  (a,b)). The scan makes the spiking region visually explicit. *Pitfall locked in:*
  the naive (a=0.7, b=0.8) choice does **not** oscillate — the notebook must use
  verified spiking parameters and show the scan that justifies them.

### NB §10 — The relaxation limit: the flagship (§7.4, §9.1)
- Markdown: the full explicit computation of §7.4 — slow manifold
  `y_s = x/[μ(1−x²)]`, slow-flow time `t_R = μ(3/2 − ln2)`, jump matching
  `C = 2μ/3`, the root `x=2` of `(x−2)(x+1)²=0`, and the boxed result
  `T(μ) = μ(3 − 2 ln 2) + O(1)`, amplitude → 2.
- Code (the centerpiece numeric):
  1. **Regime table.** For μ ∈ {0.05, 0.1, 1, 10, 50, 100}: Radau (rtol=atol=1e-10),
     measure T by zero-crossings, amplitude; build a pandas DataFrame with columns
     `T, T/μ, T/(μ(3−2ln2)), amplitude`; compare to the essay's table.
  2. **Stiffness contrast (the §9.1 payoff).** At μ=100: (a) Radau integrates one
     period in O(μ) steps; (b) an explicit RK4 at the stability-limited
     `dt = 2.785/(3μ)` **diverges in the |x|<1 strip** (print the step index where
     it blows up). Plot both: Radau's smooth cycle vs RK4's runaway.
  3. **Asymptotic overlay.** Plot measured `T/μ` vs `1/μ` for the table, with the
     horizontal line `3−2ln2 = 1.61371`; show convergence from above.
- **Expected (measured, matches essay table to 4–5 digits):**
  T(10)=19.078, T(50)=82.508, T(100)=162.837; T/μ = 1.908/1.650/1.628 → 1.6137;
  amplitude 2.014/2.003/2.001; RK4 diverges at step 4970 (t≈46).
- This single section is the reason the notebook exists: it ties the explicit
  asymptotics of §7.4 to the O(μ²)-vs-O(μ) stiffness of §9.1 in one figure.

### NB §11 — Stochastic Liénard: Boltzmann density & Kramers escape (§8)
- Code (SDE, Euler–Maruyama, `QUICK`-scaled):
  1. **Stationary density.** Linear friction f≡γ, double-well V; run the SDE
     `dx=v, dv=−(γv+V′)dt+σ dW`; discard burn-in; histogram the x-marginal;
     overlay the Boltzmann `∝ exp(−V/D)`, `D=σ²/(2γ)`; report the correlation.
  2. **Kramers escape.** Double-well V (ΔV=0.25); for D ∈ {0.15, 0.25, 0.40, 0.60}
     count barrier crossings with a **hysteresis well-state counter** (state flips
     L→R only when x crosses 0.2 from the left); report rate vs D; plot
     `ln(rate)` vs `−ΔV/D` and fit the slope (theory: 1).
- Markdown: state the underdamped caveat — the x-marginal is **not** exactly
  Boltzmann when inertia is present (corr ≈ 0.88, not 1.0); the overdamped limit
  `γ→∞` is where Boltzmann is exact and the limit cycle disappears (1-D gradient
  flow ⇒ no isolated periodic orbit).
- **Expected (measured):** corr(empirical, Boltzmann)≈0.88; Kramers rate rises
  monotonically with D (0.0167→0.0317→0.0492→0.0692), `ln rate ≈ −ΔV/D + const`,
  fitted slope ≈ 1.1. *Pitfalls locked in:* (a) a naive "crossed the barrier"
  detector that checks `x<−0.5 AND x>0.5` on consecutive samples can never fire
  (x moves continuously) — must use a hysteresis state counter; (b) short runs bias
  to one well — must discard burn-in and run long enough to visit both wells.

### NB §12 — Structure-preserving discretization (§9.4)
- Markdown: a Liénard system is not Hamiltonian (`div=−f≠0`), so strict symplectic
  integration doesn't apply; but the **energy identity (4)** still prescribes a
  splitting — conservative `(x,v)→(v,−V′)` + dissipative `(0,−f v)`.
- Code: conservative Duffing (δ=0, V=−x²/2+x⁴/4); long run (T=1000, h=0.05, N=20000):
  **leapfrog (symplectic)** vs **explicit RK4**; plot `E(t)` for both; report the
  final |ΔE|.
- **Expected (measured):** leapfrog |ΔE|=1.5e−7 vs RK4 8e−6 — leapfrog conserves
  ~50× better and stays bounded (bounded energy error) while RK4 drifts.
- Markdown: name the Biryukov benchmark (NB §7) as the test bed; cite Hairer–Lubich–
  Wanner [29].

### NB §13 — What this notebook does *not* compute (§10, §6.3)
- Markdown only: the open problems (n=5 cyclicity, exact H(n), cyclicity of
  graphics / alien derivations, canard-explosion uniform theory, certified numerics
  in the stiff regime, non-perturbative stochastic bounds). Frame each as "the
  numeric hook the notebook *would* use if it were treatable" (e.g., certified
  numerics → build on the reduced slow-flow of NB §10, not brute-force IVP).

---

## 3. Generation & execution mechanics

1. **Source of truth = builder script.** Write `build_lienard_notebook.py` using
   `nbformat.v4` (`new_notebook`, `new_markdown_cell`, `new_code_cell`). Each
   notebook section = a small function returning a list of cells; the script
   concatenates them and `nbformat.write` to
   `LienardPotential-main/LienardPotential/lienard_potential_notebook_v0.ipynb`.
   This mirrors your part-file-then-concatenate convention and keeps the notebook
   diffable/versionable.
2. **Execute in place** to bake in outputs and prove it runs:
   ```
   ./.venv/Scripts/python.exe -m jupyter nbconvert \
     --to notebook --execute --inplace \
     --ExecutePreprocessor.timeout=900 \
     LienardPotential-main/LienardPotential/lienard_potential_notebook_v0.ipynb
   ```
   (`matplotlib` backend set to `Agg` in NB §0 so headless execution works.)
3. **Open in JupyterLab** (already installed):
   ```
   ./.venv/Scripts/python.exe -m jupyter lab \
     LienardPotential-main/LienardPotential/lienard_potential_notebook_v0.ipynb
   ```

**File layout (all under the essay's directory):**
```
LienardPotential-main/LienardPotential/
  lienard_potential_essay_v3.md          (input, untouched)
  build_lienard_notebook.py              (builder script — source of truth)
  lienard_potential_notebook_v0.ipynb    (generated + executed)
```

---

## 4. Acceptance criteria (verification checklist)

The notebook is **done** when, on a clean `--execute` run, every one of the
following holds (each is a printed assertion/line in the notebook, not a hope):

- [ ] **A1** Sympy: Liénard transform and energy identity both simplify to the
      claimed forms (`True`).
- [ ] **A2** Averaging: W₂, W₄, W₆ = π, 3π/4, 5π/8; I(h) zero at h=2 ⇒ radius 2.
- [ ] **A3** Return map (NB §4): exactly one zero of the displacement; multiplier<1.
- [ ] **A4** Biryukov (NB §7): three ICs converge to one cycle, amp≈2.5012, T≈6.5658.
- [ ] **A5** Duffing (NB §8): equilibria 0,±1; V(±1)=−0.25; separatrix E=0.
- [ ] **A6** FitzHugh–Nagumo (NB §9): spiking case period≈221, spike≈4; the scan
      correctly marks (0.7, 0.8) as *non*-spiking.
- [ ] **A7** **Flagship (NB §10):** T(10)=19.078, T(50)=82.508, T(100)=162.84
      (±0.01); T/μ→1.6137 from above; amplitude→2; RK4 diverges at μ=100.
- [ ] **A8** Stochastic (NB §11): Boltzmann corr≈0.88; Kramers rate monotone in D,
      ln-rate slope≈1.
- [ ] **A9** Structure-preserving (NB §12): leapfrog |ΔE| < RK4 |ΔE| by ≥10×.
- [ ] **A10** Notebook executes end-to-end with no errors in ≤ ~3 min (or ≤ ~30 s
      with `QUICK=True`).

---

## 5. Risks & mitigations (all discovered during planning)

| Risk | Mitigation (baked into the plan) |
|---|---|
| **Half-period bug** — measuring T between wrong crossings gives T/μ≈0.81 instead of 1.61 | `period_from_crossings` enforces same-direction crossings on the settled window; A7 pins the correct values. |
| **FH-N non-spiking params** — (a=0.7, b=0.8) settles to equilibrium, no oscillation | Use verified spiking params (a=0, b=0); include the parameter *scan* that makes the spiking region explicit; A6 asserts (0.7,0.8) is non-spiking. |
| **RK4 instability** — explicit methods blow up in the \|x\|<1 strip at large μ (this is the point, but must not crash the notebook) | Wrap the RK4 demo to *detect and report* divergence (finite-check) rather than let it raise; it's a demonstration cell. |
| **Kramers detector** — "x<−0.5 AND x>0.5 on consecutive samples" never fires (continuous x) | Hysteresis well-state counter (state flips only on crossing 0.2 from the prior well); A8 asserts monotonic rate. |
| **Stochastic burn-in / one-well bias** — short runs sample only one well | Discard first half; run long enough to visit both wells; A8 reports both-well coverage. |
| **Boltzmann overclaim** — underdamped x-marginal ≠ exact Boltzmann | State the caveat explicitly; report corr≈0.88, not 1.0; reserve "exact Boltzmann" for the overdamped limit. |
| **Runtime blowup** — 1.5M-step SDE + 4-run Kramers sweep | `QUICK` flag scales step counts/dt; `--ExecutePreprocessor.timeout=900`; A10 caps total runtime. |

---

## 6. Immediate next step (on your go-ahead)

1. Write `build_lienard_notebook.py` (builder, ~13 sections as above).
2. Generate `lienard_potential_notebook_v0.ipynb`.
3. Execute it in place with `jupyter nbconvert --execute`; confirm A1–A10.
4. Hand off the executed notebook for JupyterLab.

*All headline numbers in this plan were produced by real runs in the workspace
`.venv` during planning — none are assumed.*
