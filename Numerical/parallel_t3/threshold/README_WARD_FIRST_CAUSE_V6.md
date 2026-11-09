# Stream 3 v6: isolate the first non-invariant gauge tensor

Read-only diagnostic against the **source-backed saved** rough Stream 1/2 inputs. All original data and source files are left untouched. This is not a matching theory, replacement CLEBSCH map, or physical completeness proof.

## Evidence motivating the check

- Saved `C_psi2phi2_tree_mf` has **zero LLHH** for all 16 models; this says nothing about whether its **LLSS** sector satisfies SU(2).
- At `M_S`, 14 LLHH sectors are holomorphic HH (correct U(1) charge) but violate SU(2). B(-1)/C(-1) also have non-holomorphic sectors and fail U(1).
- In the inspected 2026-11-10 GitHub source, `Numerical/parallel_t3/rough_inputs/builder.py:quartic_unit_tensors` generates a conditional `lambdaT3` for `(2,2)`, `(3,1)`, `(3,3)` using an index contraction that is **not SU(2) invariant** in `wilson.basis.build_basis` real coordinates. Direct complex-representation contraction gives relative norms ~1.0, ~1.632993, ~1.632993, while `(1,3)` is invariant. This source issue alone does not establish the cause in A(-4), nor the B(-1)/C(-1) hypercharge violation.

## One diagnostic command (PowerShell)

From the project root, after extracting the v6 zip:

```powershell
.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_rough_first_cause_v6.py
.\.venv\Scripts\python.exe scripts\AuditT3RoughFirstCauseV6.py --models A_-4 --with-beta
.\.venv\Scripts\python.exe scripts\AuditT3RoughFirstCauseV6.py
Import-Csv output\parallel_t3\rough_results\rough_fermion_first_v1\_diagnostics\first_cause_v6\first_cause_summary.csv |
  Format-Table model,MF_Wilson_SU2,MF_Yukawa_SU2,MF_quartic_SU2,MF_saved_beta_SU2,MS_Wilson_SU2,bad_source_quartic_unit_names -AutoSize
```

The script verifies the Stream 2 input SHA-256 recorded for the Stream 1 canonical `background.npz` against the **current** file. If the background changed, it blocks rather than mixing old Wilson and new running backgrounds. It validates basis labels, complete tensor dimensions and both physical scales.

Optional `--with-beta` calls the *existing* `Numerical.parallel_t3.wilson.adm.beta` once at the saved `M_F` endpoint (expensive for the 45-Weyl tensor) and reports each contribution from the seven implemented ADM channels. No UV, threshold or Wilson ODE is re-run. The default scan uses the already-saved `beta_C_at_MF_16pi2` and `beta_C_at_MS_16pi2` instead.

## Output

`output/parallel_t3/rough_results/rough_fermion_first_v1/_diagnostics/first_cause_v6/`

- `<slug>_first_cause.json`: input fingerprints; full MF/MS Wilson/Yukawa/quartic SU(2), U(1), SU(3) Ward norms; optional saved-vs-recomputed beta disagreement; source-coded quartic **component-unit** Ward tests; existing MF/MS LLHH charge sectors.
- `first_cause_summary.csv` and `first_cause_summary.json`: flattened all-16 diagnostics and specific first non-invariant sectors.

Interpretation: if `MF_Wilson_SU2` is large, the assumed formal Matchete LLSS → canonical Weyl/real-field map must be corrected **before any ADM conclusion**. If `MF_Wilson_SU2` is tiny but `MF_quartic_SU2` is large, fix the native quartic CG/dictionary first. If both inputs and Yukawa are invariant but `MF_saved_beta_SU2` is large, inspect the full candidate ADM and its channel breakdown. Even if all checks pass, the result is still an unverified rough numerical benchmark. Never overwrite an old source with an unproved projected tensor or call the projected result a gauge-invariant SMEFT Wilson coefficient.
