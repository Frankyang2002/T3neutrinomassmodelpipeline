# Track B — Independent direct UV A/D/E nine-model batch

## Source boundary and inventory

GitHub read-only repository: `Frankyang2002/T3neutrinomassmodelpipeline`, default branch `main`; inspected recursive tree at commit `d2163e6bb95109c0fac419f7e41f335f0d347b62` (GitHub timestamp 2026-11-03T01:48:36Z). **The actual Windows local HEAD and uncommitted state have not been observed.** The delivery does not alter GitHub or shared source. The only locally accessible project materials here are the uploaded reference PDFs; no Windows repository or executable Matchete kernel is available in this execution container. This is a source-inspected, locally executable handoff; no unattended certification is claimed.

## Ordered work stages (fixed)

1. **B0 — Provenance and all-nine inventory.** Verify current source versions, restricted model representations, available native witnesses and benchmark inputs. Establish a nine-row evidence ledger. **Status: DONE for GitHub source/tree/benchmark and nine-row inventory; local Windows HEAD verification pending.**
2. **B1 — Native source normalisation.** Extract or independently verify every model's charge-selected Yukawa CGs, scalar mixing CGs, full portal invariant decomposition, relevant scalar mass blocks and Dirac fermion mass conventions. **Status: COMPLETE conditionally using archived actual Matchete CG registries, the full invariant identity span, all charge channels and source Hessians; fresh nine-model live scalar/Yukawa/Dirac-mass actions PENDING.**
3. **B2 — Independent nine-model UV calculation.** Derive the charge-conserving contributions and compute all nine neutrino matrices without importing the provisional UV loop engine. **Status: NINE archived-native-CG zero-portal and separate conditional benchmark-identity matrices evaluated; local nine-model Python recalculation implemented, fresh-native run PENDING.**
4. **B3 — Numerical validation.** Test UV pole cancellation, mass positivity, matrix symmetry, Takagi singular values, complex rephasing covariance, independent loop-integral agreement and comparison with existing frozen matrices only AFTER independent calculation. **Status: source-mass positivity, pole unitarity, matrix symmetry, Takagi values, adaptive quadrature vs 75-digit B0 and comparison tested; synthetic three-family complex/charged engine tests PASS; nine real-benchmark complex/native reruns PENDING.**
5. **B4 — Single consolidated deliverable.** Package full source files, nine result reports, model-specific certification status, one batch summary, one ZIP and exact PowerShell commands. **Status: DONE: one ZIP with source, nine reports, consolidated batch summary, timestamped evidence and exact PowerShell commands. Physical certification remains PENDING for all nine.**

## Nine-model evidence ledger (B0→B4 snapshot)

| Model | Restricted `(dS1,dS2,dF,alpha)` | Native source | Scalar/fermion normalisation | Independent matrix | Takagi descending (eV) | Frozen residual | Eligible |
|---|---|---|---|---|---|---|---|
| A_-4 | (1,3,2,-4) | archived actual CG, live action pending | portal identity conditional; indexed Dirac mass/spinor phase pending | two source points evaluated | 0.043600644, 0.0075827281, 0.00047448293 | 8.87e-12 | no |
| A_-2 | (1,3,2,-2) | archived actual CG, live action pending | portal identity conditional; indexed Dirac mass/spinor phase pending | two source points evaluated | 0.048560209, 0.011008458, 0.0072147535 | 8.53e-12 | no |
| D_-2 | (3,1,2,-2) | archived actual CG, live action pending | portal identity conditional; indexed Dirac mass/spinor phase pending | two source points evaluated | 0.048482372, 0.010997754, 0.0072141749 | 7.54e-12 | no |
| D_+0 | (3,1,2,0) | archived actual CG, live action pending | portal identity conditional; indexed Dirac mass/spinor phase pending | two source points evaluated | 0.048549903, 0.011007039, 0.0072146703 | 7.29e-12 | no |
| D_+2 | (3,1,2,2) | archived actual CG, live action pending | portal identity conditional; indexed Dirac mass/spinor phase pending | two source points evaluated | 0.04359482, 0.0075817378, 0.00047473317 | 7.84e-12 | no |
| E_-4 | (3,3,2,-4) | archived actual CG, live action pending | portal identity conditional; indexed Dirac mass/spinor phase pending | two source points evaluated | 0.038488455, 0.0074169371, 0.0033193971 | 1.5e-11 | no |
| E_-2 | (3,3,2,-2) | archived actual CG, live action pending | portal identity conditional; indexed Dirac mass/spinor phase pending | two source points evaluated | 0.042469478, 0.0086788288, 0.0047141911 | 2.12e-11 | no |
| E_+0 | (3,3,2,0) | archived actual CG, live action pending | portal identity conditional; indexed Dirac mass/spinor phase pending | two source points evaluated | 0.042515357, 0.0086485894, 0.0046435843 | 2.58e-11 | no |
| E_+2 | (3,3,2,2) | archived actual CG, live action pending | portal identity conditional; indexed Dirac mass/spinor phase pending | two source points evaluated | 0.03863356, 0.0074424848, 0.0033259668 | 2.03e-11 | no |

No A_+0 mass or normalisation is imported. The published GitHub artifacts for A_+0/B_-1/C_-1 are not native witnesses for these nine models.

## Acceptance gates

For each model: exact source hashes and restricted representation verified; native action/yukawa and conjugate orientation measured; complete SU(2) invariant decomposition and charge-selected components measured; physical Hessian positive; Dirac mass and spinor conventions verified; a source-only self-energy assembled and finite; independent numerical loop representation agrees; complex basis covariance and Takagi checked. A small residual to frozen results is *not* sufficient for certification. Conditional identity portal and zero-portal results must be distinctly labelled. A common 100/110/120 TeV mass list is used only when checked against indexed native mass terms and source benchmark. Any missing gate is `PENDING`, not zero, not a fitted factor.


## Consolidated B4 result and unresolved gates

The nine matrices in `ARCHIVED_NATIVE_NINE_MODEL_SNAPSHOT.json` and `reports/` were evaluated directly from archived native Matchete CG components and source quadratic/mixing conventions, **not** fresh native-action-certified physical matrices. The local Python runner reproduces these and extends to real-source complex benchmark tests on Windows. A new timestamped `output/general_observable/physical_normalisation/parallel_B_ADE_.../` directory in this single ZIP stores full independent snapshot reports; later local executions must make fresh directories rather than overwriting this evidence.

**Exact failing physical gates for every row:** no fresh model-specific assembled `CheckLagrangian` + scalar quadratic + flavour-indexed Dirac `MF[r]` witness yet; unresolved installed numerical Matchete `GammaCC`/Weyl relative convention; and physical identity-portal RGBeta-to-native dictionary not established. A frozen residual below `3e-11` cannot resolve these gates. Local `git status` and source SHA checks precede acceptance of a new native witness.

## B4 runtime status update (generated this run, not a source edit)

Run directory: `output\general_observable\physical_normalisation\parallel_B_ADE_20261103_133650_813175`. B0 source inventory checked; B1 and B2 numerical gates per row below; B3 numerical controls tested for completed rows; B4 consolidated reports written. Physical eligibility remains false until unresolved native source and spinor gates are independently measured.

| Model | Native action | Matrix | Takagi descending (eV) | Comparison residual | Eligible |
|---|---|---|---|---|---|
| A_-4 | pending | blocked: A_-4: invalid fresh native action witness | pending | pending | no |
| A_-2 | pending | blocked: A_-2: invalid fresh native action witness | pending | pending | no |
| D_-2 | pending | blocked: D_-2: invalid fresh native action witness | pending | pending | no |
| D_+0 | pending | blocked: D_+0: invalid fresh native action witness | pending | pending | no |
| D_+2 | pending | blocked: D_+2: invalid fresh native action witness | pending | pending | no |
| E_-4 | pending | blocked: E_-4: invalid fresh native action witness | pending | pending | no |
| E_-2 | pending | blocked: E_-2: invalid fresh native action witness | pending | pending | no |
| E_+0 | pending | blocked: E_+0: invalid fresh native action witness | pending | pending | no |
| E_+2 | pending | blocked: E_+2: invalid fresh native action witness | pending | pending | no |
