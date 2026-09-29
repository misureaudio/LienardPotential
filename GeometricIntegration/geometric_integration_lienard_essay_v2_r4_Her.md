The builder just executes the notebook — I need to read the notebook itself to compare the numerical methods.I read the builder (it just executes + validates the notebook) and the notebook itself, then cross-checked its numerical methods against the attached v2 scripts, essay, and remediation plan.

## Verdict

Largely **consistent** — the v2 remediation implements exactly what the notebook diagnosed (correct `c = -h/2·J₂₁`, fail-loud guard, script-generated Markdown tables, bounded-but-distorted IM narrative). The Jacobian sign, Cramer's update lines, stability function `R(z)=(1+z/2)/(1−z/2)`, `λ_max=3μ`, the Radau reference (`rtol=1e-10`, `max_step=0.2`), the Strang split, and the IM/Strang/Radau results all agree between the notebook and the attached files.

But there are a few **numerical-method inconsistencies** in the verification path worth flagging:

### 1. Newton iteration cap (`itmax`) — the one real conflict
- Notebook "correct" Newton (cell 8) and `computation3.py`: **`itmax=200`**.
- `computation4.py`: **`itmax=2000`** (`def imid_step(..., itmax=2000, guard=1.0)`).

`computation4.py` is the script that generates the §6.3 essay table, yet its corrected IM uses a 10× larger Newton cap than the notebook's *verified* "correct IM." Since the jump state is precisely where convergence is marginal, this can change whether a stiff step resolves. This is the likely source of the small `max|x|` drift at `h=0.5`: notebook conclusion says **3.0**, the remediation plan says **3.03**, the v2 essay table says **3.05** — three different values for the same 3-period run.

### 2. Fail-loud guard: in the deliverable, absent from the verified method
- v2 `computation3.py`/`4.py`: `raise RuntimeError` if `res > guard=1.0`.
- Notebook "correct" Newton (cell 8): **no guard** — it silently returns the last iterate (`conv=False`), and `run_orbit` (cell 13) discards the flag via `[:2]`.

So the notebook proves "bounded-but-distorted" with a *guard-less* IM, while the shipped scripts use a *guard-ed* IM. At any state where the correct Newton fails to converge with residual > 1.0, the notebook's IM would continue silently (→ "jump unresolved") but the v2 scripts would **abort** (→ `UNRESOLVED`). They only agree because the residual happens to stay < 1.0 at the jump — a narrow margin the verification doesn't actually check.

### 3. Independent root-solver cross-check: in the notebook, not in the v2 scripts
The notebook's Method II (cell 8) validates the corrected Newton against `scipy.optimize.root` (hybr) at multiple stiff states. The v2 `computation3/4.py` drop this entirely — they rely solely on the hand-coded Newton + residual guard. So the notebook's key *independent* verification of the corrected Jacobian is not replicated in the deliverable (the remediation plan offered "fall back to `scipy.optimize.root`" as Option B, but the scripts took Option A and kept no cross-check).

### By-design differences (not inconsistencies)
- The notebook's `imid_b` (buggy, `tol=1e-13, itmax=100`) reproduces the **v1** scripts, not the attached v2 ones — expected, it's the forensic target.
- v2 tightened `tol` from `1e-13`→`1e-14` per the plan — intentional.

**Bottom line:** the physics narrative and the corrected Jacobian are consistent; the friction is all in the Newton *solver configuration* — `itmax` (200 vs 2000), the presence/absence of the fail-loud guard, and the missing independent `scipy.optimize.root` check. If you want the essay's §6.3 numbers to be exactly reproducible by the notebook's verification, align `computation4.py` to `itmax=200` (or re-verify at 2000) and add the guard behavior to the notebook's "correct" IM so both fail identically.