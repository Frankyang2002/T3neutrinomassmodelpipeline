# Batch 2: repair the A(-4)/D(-2) raw Laurent audit

Date: 2026-11-10. Source inspected **read-only**: `Frankyang2002/T3neutrinomassmodelpipeline`, GitHub `main` commit `ab8f5452a0887f155d063ac12a790c6c9c1134be`.

## Verified failure

The uploaded GitHub outputs in
`output/parallel_t3/full_general_batch2/thesis_full_general_batch2_laurent_v1/{A_-4,D_-2}`
show `Blocked_LoopBookkeeping` with reason `Correction is not exactly hbar-linear and nonzero` for both models. The exact native raw one-loop export paths under `thesis_full_general_batch2_raw_v1` contain:

| Model | Top-level additive InputForm terms | `hbar` terms | no-`hbar` terms |
|---|---:|---:|---:|
| A_-4 | 135 | 131 | 4 |
| D_-2 | 1,200 | 1,194 | 6 |

The 4 and 6 zero-loop-order terms are pairs of scalar quadratic mass-coordinate terms. They use different native Matchete mass couplings (`CNewScalar12` / `MS1^2`, `CNewScalar22` / `MS2^2`, and in D(-2) the `CT3Scalar2{Minus,Plus}2` / `T32{Minus,Plus}Mass^2` pairs). They are *not* automatically identified or set to zero by this repair; whether they represent a valid tree-level renormalisation or a source mismatch requires checking separately.

The raw InputForm also uses **`ϵ`** (`\[CurlyEpsilon]`) rather than only `ε` (`\[Epsilon]`). The old auditor did not match this regulator spelling.

## Replacement files (full files, not patches)

- `Lagrangian/RunT3Batch2LaurentAudit.wl`: full replacement. Separates the full source exactly into its `hbar^0` and `hbar^1` pieces, rejects `hbar^2` or unrelated tree-level operators, separately exports the unreconciled mass remainder and hbar-linear coefficient as WXF, recognizes both Greek epsilon variants, and then performs the formal Laurent audit. Reuses source loops if no unevaluated named loop heads remain. **No coefficient is silently discarded.**
- `Numerical/parallel_t3/wilson/general_batch2_laurent.py`: full replacement, validates both preserved components, term counts, new v2 schema, operator/checksum coverage and the absence of physical matching claims.
- `scripts/RunT3ThesisGeneralBatch2Laurent.py`: full replacement, surfaces the native `Blocked_*` Wolfram stage rather than hiding it behind `Blocked_RuntimeError`.
- `tests/test_t3_thesis_general_batch2_laurent.py`: full replacement; 17 synthetic tests, no Wolfram or physics certification.
- `docs/T3_BATCH2_LAURENT_REPAIR_2026_11_10.md`: this handover.

No changes to Batch 1, source raw Matchete matching, physics constants, or pre-existing output directories.

## Windows execution

From the project root `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
Set-Location $repo
$zip = 'C:\path\to\T3_FullGeneral_Batch2_Laurent_Fix_2026_11_10.zip'
Expand-Archive -LiteralPath $zip -DestinationPath $repo -Force

python -X utf8 -m unittest discover -s tests -p test_t3_thesis_general_batch2_laurent.py -v
python -X utf8 -m scripts.RunT3ThesisGeneralBatch2Laurent --all --preflight
python -X utf8 -m scripts.RunT3ThesisGeneralBatch2Laurent `
  --all `
  --run-id thesis_full_general_batch2_laurent_v2 `
  --accept-formal-audit-only
```

**Use the new `v2` run ID.** The previous v1 reports remain untouched. Do not rerun either RG segment or the raw Matchete F-only first threshold.

## New outputs and interpretation

`output/parallel_t3/full_general_batch2/thesis_full_general_batch2_laurent_v2/<model>/`

- `unreconciled_hbar0_mass_remainder.wxf`: **unresolved** tree-order scalar mass-coordinate difference. Its presence is recorded, not silently cancelled.
- `raw_hbar1_coefficient.wxf`: complete formal order-one-loop factor-stripped coefficient, not numerical.
- `evaluated_full_hbar1_coefficient.wxf`: source-preserving formal result after native loop-function processing when needed.
- `explicit_pole_residue_unsubtracted.wxf`, `epsilon_zero_regular_unsubtracted.wxf`: produced only when explicit Laurent extraction is successful.
- `laurent_audit.json`: per-model native stage, regulator, term counts, exactness flags and SHA256 of WXF outputs.
- `manifest.json` and `wolfram_laurent.log`: Python validation outcome and Wolfram log.

On success, the audit status is `LoopFunctionsEvaluated_ExplicitLaurentSeparated_TreeRemainderPreserved_AuditOnly`, and Python's status is `FullOneLoopFormalLaurentAudited_TreeMassRemainderPreserved_NotFiniteWilsons`.

**Not implemented by this repair:** UV-hard minus EFT1-soft subtraction, local MSbar counterterms and pole cancellation, the physically correct Matchete-to-RGBeta and real-scalar basis conversion, finite numerical `C_ijab(MF)`, EFT1 Wilson RGEs, scalar matching, SMEFT, masses or chi-square. The regular part of the Laurent expansion is NOT a renormalized Wilson coefficient. D(-2)'s extra self-bilinear and all five additional quartics are preserved by the original numerical boundary inputs.

## Verification performed for this package

- Current GitHub source, both failing `laurent_audit.json`, both `wolfram_laurent.log`, and exact raw InputForm exports inspected.
- Matched operator counts and hbar-zero mass terms verified from raw source.
- Native Wolfram algebra checks on small representative A/D mass terms, hbar polynomial splitting and epsilon symbol recognition passed in an external Wolfram kernel.
- Python compilation and **17/17 offline unit tests passed**.

**Not done:** The full A/D local Matchete Laurent execution was not run in this environment; subsequent failure at loop-function or Laurent steps must be diagnosed from the new v2 reports, not reclassified as success.
