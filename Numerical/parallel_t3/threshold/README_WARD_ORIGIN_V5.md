# Workstream C v5 — trace the source LLHH Ward obstruction

Scope: **source-only diagnosis** of the existing Stream 2 `C_psi2phi2_ms` tensor.
The original simultaneous-scalar-matching algorithm, ADM implementation,
SMEFT runner, and all input files remain untouched.  No tensor is silently
symmetrised, repaired, reindexed, or truncated.  There is **no second EFT**.

## Evidence motivating this tool

The local Windows v4 regression suite passed 131 tests, but the genuine
A/B all-16 run stopped at the existing C preflight: LLHH Weinberg projection
residuals 0.9479 to 1.0.  In 14 models the U(1)_Y Ward residual was near
machine precision while all three SU(2) Ward residuals were generally order
one.  The two mixed-doublet points B_-1 and C_-1 additionally had U(1)_Y Ward
residuals 1.3987 and 1.5209.  Those are **source failures in the declared
basis**, not a scalar-bubble precision tolerance issue.

## Physical algebra of the diagnostic

With the *declared* Higgs coordinates `H=(R+iI)/sqrt(2)` and ordinary SM
lepton doublet hypercharge -1/2, convert both Higgs legs of the full LLHH
coefficient into complex components `H1,H2,H1*,H2*`.  The conversion matrix
is unitary (2 real Higgs components per complex field), so the amplitude norm
splits exactly as

```
||C_LLHH||_F^2 = ||C_LLHH[HH]||_F^2
               + ||C_LLHH[H H* / H* H]||_F^2
               + ||C_LLHH[H* H*]||_F^2.
```

Their total operator hypercharges (including `LL`) are respectively 0, -1,
-2, giving an independent identity:

```
||delta_Y C_LLHH|| / ||C_LLHH||
    = sqrt( ||C[H H* / H* H]||^2 + 4 ||C[H* H*]||^2 ) / ||C_LLHH||.
```

The code verifies this algebraic identity independently against the direct
U(1)_Y Ward evaluation.  For the **allowed HH sector only**, the code tests
all three SU(2) generators and the sole Weinberg singlet, never using a
convenient neutral-component value as a substitute.  No new operator or
finite matching contribution is invented.

`C_psi2phi2_tree_mf` (if present) is checked separately.  A zero order-zero
LLHH at MF is explicitly *uninformative* about the running origin; it does
not prove a source-derived anomalous-dimension bug.  A nonzero pure MF tree
and impure MS endpoint localises the problem **after the saved tree endpoint**,
which could include hard MF terms, Wilson initialisation, its input trajectory,
its RGE or a basis map. It is **not** automatically proof that the ADM is wrong.

## Installation / tests / origin readout

Run from the repository root:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Join-Path $env:USERPROFILE 'Downloads\t3_rough_stream3_ward_origin_v5.zip'
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_rough_ward_origin.py

# Consumes the `_adapted_inputs` already written by the successful v4 adapter.
# Does NOT rerun UV, MF, EFT1, the MS scalar threshold, or SMEFT RGEs.
.\.venv\Scripts\python.exe scripts\AuditT3RoughWardOrigin.py --models A_-4
.\.venv\Scripts\python.exe scripts\AuditT3RoughWardOrigin.py

Import-Csv output\parallel_t3\rough_results\rough_fermion_first_v1\_diagnostics\origin_v5\origin_summary.csv |
  Format-Table model,MS_charge_category,MS_fraction_HH,MS_fraction_H_Hbar,MS_fraction_Hbar_Hbar,MS_HH_Singlet_residual,MF_orderzero_charge_category -AutoSize
```

Per model `output/parallel_t3/rough_results/rough_fermion_first_v1/_diagnostics/origin_v5/<slug>_origin.json` is written; `origin_summary.csv` and `origin_summary.json` aggregate all models.  These are C-owned diagnostic files only.  A missing or invalid endpoint is recorded as a blocker instead of being fabricated.

## An explicitly projected numerical calculation (separate output root)

The previous v4 matching runner already supports an explicitly labelled
`NonWeinbergProjectedDiagnostic` mode.  This is **not** a valid gauge-invariant
SMEFT mass prediction.  If a rough numerical engineering smoke test is useful:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3RoughAll16.py `
    --models A_-4 --allow-projection-residual --overwrite `
    --output output\parallel_t3\rough_results\projection_diagnostic_v1

# Only after checking the A(-4) result, process all 16 diagnostic projections.
.\.venv\Scripts\python.exe scripts\RunT3RoughAll16.py `
    --allow-projection-residual --overwrite `
    --output output\parallel_t3\rough_results\projection_diagnostic_v1
```

The strict required output root `rough_fermion_first_v1` is left untouched by
this diagnostic experiment.  Neither this script nor the opt-in diagnostic
mode is a replacement for repairing the source EFT1 tensor or basis bridge.

## Checks and limitations

- Reuses source C `extract_lepton_sectors` and `master_weinberg_template`,
  and the existing checksum-verified A endpoint loader.
- The source's full tensor and separate C5 are never rewritten.
- Algebraically classifies *allowed-hypercharge but non-SU2* versus
  *wrong-hypercharge* LLHH contributions.
- Does not claim the underlying signed GammaCC/Weyl normalisation, full
  ADM, source/native quartic dictionary or complete one-loop hard matching.
- Standalone package test: **10/10 passed**; full Windows regression suite
  and user model outcomes must be verified in the user's local checkout.
