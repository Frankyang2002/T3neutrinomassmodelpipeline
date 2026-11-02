# B2 batch 3 — handoff, exact literal-index Ward witness

## Source and scope

- GitHub `main` inspected: `e2a5def5326cc452f668ef593bd81ae13fdff30b`.
- RGBeta model source Git blob: `7af7a4ddf99a7e38210543a21b3b1cda203eaf02`; native formal file Git blob: `34bf2f1d7b325cab0d5162b235ea2df0e8034c78`.
- `analytic_reference.json` is **pure analytic chosen-convention evidence**, not a Windows-native execution report. `RunPhysicalB2RGBetaWardWitness.py` creates a new, locally source-guarded report. It accepts a different current HEAD if and only if all six previously reviewed source blobs remain identical and clean.

## Measured and calculated

- Exact SU(2) Ward residual norm squares of literal native-formal B/C `delta*delta` scalar mixing: `(3,0,3)`, **non-invariant under declared literal fund/bar conventions**.
- Exact residual for a chosen symmetric `delta*epsilon` mixing candidate: `(0,0,0)`.
- For B-class `y1`, the formal `eps` on `(Bar[S1],l)` has `(2,0,2)`; a chosen mixed-index `delta` has `(0,0,0)`.
- For B-class `y2`, `eps` on `(S2,l)` has `(0,0,0)`.
- In B/C for each alpha `-3,-1,+1`, the literal mixing candidate violates electric-charge pairing after `H0H0`; the chosen invariant does not. Six analytic checks, five tensor rows.
- Tests: **21 exact/synthetic pytest tests passed in the authored container**, including Git head drift acceptance and real source drift rejection. **No Wolfram executed here**.

## Interpretation and integration

The physical meaning of RGBeta's `Bar`, `del`, `eps`, `delS2` index semantics is **unverified**. The package may apply a pseudoreality/dual-index map. Therefore this is a **conditional warning**, not a proof that RGBeta's actual action is wrong. No source patch is authorised. Any Matchete tensor ratio to the chosen Ward candidate is **not** a physical coupling map.

- `physical_source_certified = false`, **0/16**.
- Physical `y1,y2,MF,lambdaT3` cross-package ratios = **unresolved**.
- **Safe to integrate as gauge-invariance warning / test harness:** yes, with conditional labels.
- **Safe to integrate as physical normalisation or EFT input:** **no**.

## Required next data

1. Locally run the Python witness against committed `formal_invariants.json` and old Matchete native audits.
2. Run `RunPhysicalB2RGBetaRegistryAudit.wl` in a fresh RGBeta kernel to export actual record fields and formal functions.
3. Resolve RGBeta numerical component axes + its `Bar`/epsilon index raising, then full action factors and B1 spinor phase. The symbolic warning alone cannot settle them.
