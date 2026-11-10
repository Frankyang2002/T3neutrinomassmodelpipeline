# Batch 14 — Investigate all three remaining physical EFT gaps together

## Purpose

The **numerical** A(−4), D(−2) pipeline has already run successfully through
Batch 13. Do not rerun the 18-parameter fit. This batch performs a single
source-driven investigation of the three remaining issues rather than
silently setting their unresolved coefficients to unity.

**It is a combined physics audit, not a promise that all three can be certified
from the archived inputs.** The output must remain conditional until
independent loop/operator evidence exists.

### 1. Matchete GammaCC → physical Weyl/Weinberg operator

- Validate the source `GammaCC`/`Proj[±1]` LLHH inventory:
  3 terms per chirality in A(−4), 6 per chirality in D(−2).
- Measure the **real relative multiplier** between the source-native
  Laurent-derived EFT1 logarithm and the independent full-master projection.
  Report its complex phase mismatch and residual.
- Keep `absolute_native_GammaCC_to_physical_Weyl_operator_anchored=false`.
  Algebraically, the matching-scale logarithm and pole cancellation are both
  invariant under a common operator normalization change. They cannot fix the
  absolute sign, the 1/2 operator conventions, or the physical C5 mapping.
  Required independent input: an explicit, canonically normalized spinor
  matrix-element calculation or a source-verified Matchete basis-operator
  conversion in **the same EFT/renormalization scheme**.

### 2. Complete finite one-loop matching

- Audit **every** native F-threshold ledger sector, not just the LLHH subset:
  - A(−4): 135 operators = 81 other + 48 LLSS + 6 LLHH.
  - D(−2): 1200 operators = 836 other + 352 LLSS + 12 LLHH.
- Recheck F-hard and joint-scalar-soft pole cancellation using the saved
  complex 3×3 matrices, and recomposition of the finite source-sector matrix.
- Cross-check the native scalar source counts (A: 10, D: 68) and previously
  verified pole/log and scale-partition residuals.
- Do **not** interpret these necessary checks as the complete one-loop
  MS-bar counterterm, evanescent-operator and finite operator-mixing proof.
  Full finite matching requires independently reconciled counterterms and an
  explicitly complete physical operator basis; no missing sector is set zero.

### 3. Gauge-complete EFT1 ADM and scalar invariant mapping

- Reconstruct the exact same 45-Weyl/12-real-scalar generators, SM Yukawas,
  baseline quartics and D split-real basis from the saved source states.
- Evaluate SU(3)c, SU(2)L and U(1)Y infinitesimal Ward residuals for the full
  saved F-tree tensor, full seven-channel beta, SM Yukawa tensor, and the
  scalar potential. All are independently checkable **necessary conditions**.
- Embed a test three-generation SM Weinberg operator into the full EFT1 basis,
  set all BSM potential interactions to zero, and compare the candidate full
  tensor beta against the **separate** `WeinbergRGE.sm_weinberg_beta` implementation.
  This tests several coefficients/sign conventions without assuming the BSM
  anomalous dimension is correct.
- Rebuild each archived Matchete CG scalar quartic channel at unit coupling,
  compare it to the RGBeta-potential invariant, and report the measured
  native/reference multiplier and tensor-direction residual for each supported
  alias. This is a **candidate** mapping check, not proof for missing D
  quartics or all complex phases. The two D extra HHSS tensors and the
  three scalar-only quartics still require independent basis certification.
- The published candidate seven-term anomalous dimension is not automatically
  loop-certified by Ward covariance or its SM limit.

## Installation (PowerShell)

From the repository root:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive `
  -LiteralPath "$HOME\Downloads\t3_three_physics_gates_batch14.zip" `
  -DestinationPath . -Force

.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_physical_closure_batch14.py
.\.venv\Scripts\python.exe -X utf8 -m scripts.RunT3PhysicalClosureBatch14 --all --preflight
.\.venv\Scripts\python.exe -X utf8 -m scripts.RunT3PhysicalClosureBatch14 --all --debug
```

The production pass computes gauge and native quartic probes for both models.
It reuses immutable previous output. The default output is:

`output/parallel_t3/thesis_physical_closure_batch14_v1/<model>/result.json`

`output/parallel_t3/thesis_physical_closure_batch14_v1/summary.json`

Optional `--generator-details` includes every SU(3), SU(2), U(1) residual in
JSON; the default saves their maxima. `--skip-native-quartic-map` skips the
expensive SymPy native quartic reconstruction, but then coupling-map closure
cannot be assessed. An `--output` argument selects a fresh directory.

The successful audit status will be:

`ThreePhysicsClosureGatesInvestigated_StillConditionalNotPhysicalC5`

This says the **audit completed**, not that the three physics gates passed.
If the SM-limit or Ward checks fail, look at
`gauge_Wilson_gate.max_relative_Ward_residuals` and
`gauge_Wilson_gate.SM_Weinberg_independent_limit`. If quartic normalizations
differ, check
`native_RGBeta_scalar_map_gate.archived_native_channel_comparison`.

**Never** promote `physical_C5_certified` to true based solely on these
same-source tests. Do not rerun the UV optimizer to repair a physics mismatch.

## Full replacement files (project-relative)

- `Numerical/parallel_t3/complete/two_model_physical_closure_batch14.py`
- `scripts/RunT3PhysicalClosureBatch14.py`
- `tests/test_t3_physical_closure_batch14.py`
- `README_T3_THREE_PHYSICS_GATES_BATCH14.md`

No changes to GitHub. No previous physics outputs modified.

### Minimal outstanding evidence for a physical C5

1. A normalized, chiral, signed spinor-operator crosswalk from Matchete
   `GammaCC` into `L5=(1/2) C5 (LH)(LH)+h.c.` with the project's
   `m_nu=-v² C5/2` convention.
2. A full hard-minus-EFT1-soft finite one-loop matching and counterterm
   computation in a clearly specified MS-bar and evanescent scheme, including
   operator closure and completeness.
3. A loop-level gauge-complete EFT1 anomalous-dimension derivation in the
   same Weyl/real-scalar convention, and independent Matchete ↔ RGBeta
   invariant normalization for all operators (including the five D extras).
