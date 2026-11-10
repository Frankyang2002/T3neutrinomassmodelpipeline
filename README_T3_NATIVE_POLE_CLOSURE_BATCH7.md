# T3 Batch 7 — source-preserving EFT1 counterterms and native Laurent-pole audit

## Scope and exact source revision

Prepared against `Frankyang2002/T3neutrinomassmodelpipeline` `main` at commit
`4fdae67ce2ee5f22d64284a2c83ebc4ffb3400f0` (inspected 10 November 2026).
This ZIP contains **five new files**. It does **not** modify GitHub or overwrite
Batches 1–6. It uses their saved numerical outputs and their existing trusted
intermediate-EFT, scalar matching, and SMEFT runners.

Required inputs under `output/parallel_t3/`:

- `fhard_triangle_all16_batch1_v1/`
- `strict_one_loop_all16_batch2_v1/`
- `native_finite_scenario_all16_batch4_v1/`
- `hard_soft_matching_all16_batch5_v1/`
- Optional for independent comparison: `eft1_msbar_counterterms_all16_batch6_v1/`

The current Batch 6 report has `15/16` models, and `B_-1` failed in
`write_counterterm_results` with `manifest C5: expected finite complex symmetric 3x3`.
The original Batch 6 manifest function `cm` invokes `symmetric3` on **individual
RGBeta quartic-source terms**, even though symmetry is a criterion on the final
physical C5 (and need not hold on each unsummed bookkeeping contribution).
Batch 7 retains every individual complex source matrix **without symmetrizing**,
stores its symmetric and antisymmetric parts as diagnostics, and verifies the
complete source sum against a flavour-symmetric total. This targets the observed
B_-1 failure. Whether the Windows B_-1 run now succeeds remains to be checked.

## Physics implemented

The order is fixed:

`UV -> ALL heavy fermions out at MF -> EFT1 with BOTH scalar multiplets -> BOTH scalars out at MS -> SMEFT`.

The EFT1 scalar counterterm is calculated with the already verified Batch 6
full canonical scalar tensors:

\[
\delta C_5^{\mathrm{EFT1,scalar}}\big|_{1/\bar\epsilon}
=+\frac{1}{2(16\pi^2)\bar\epsilon}
\mathcal P_5\!\left[C_{LLSS}\lambda_{HHSS}\right].
\]

The independent checks are (a) the negative of the saved bare scalar-bubble
pole, and (b) the finite B0 matching-scale derivative. All three fermions and
both scalar representations remain active at their correct stages.

For each source-resolved projected matrix `X_a`:

\[
X_a^{(S)}=\frac{X_a+X_a^T}{2},\qquad
X_a^{(A)}=\frac{X_a-X_a^T}{2},\quad
X_{\mathrm{total}}=\sum_a X_a, \quad
\left\|X_{\mathrm{total}}-X_{\mathrm{total}}^T\right\| \simeq 0.
\]

**No source-level antisymmetric tensor is deleted or modified.** The total
is rejected if it violates flavour symmetry, if any source is nonfinite,
or if its full sum does not reproduce the original total coefficient.

The native Matchete pole is reconstructed from the saved Laurent coefficient
per PL kernel and heavy fermion, grouped under the *original* native quartic
names. Batch 7 compares this raw native Laurent object with the explicitly
calculated scalar CT and with the **seven-channel candidate EFT1 anomalous-
dimension** residue. It records both raw differences, scale diagnostics, and
best-fit complex scalar normalisation ratios **for diagnosis only**. Neither
an inferred normalisation nor an invented UV counterterm enters the masses.

If the user selects the default full run, the three *conditional* C5 boundary
scenarios are independently propagated through the existing SMEFT evolution
at subtraction scales `MF/2`, `MF`, `2*MF`. At the nominal scale all three
branches must reproduce Batches 2, 4, and 5.

## What is NOT verified

1. The saved Matchete `MatchedEFT` pole has **not** been proven to equal an
   unsubtracted bare UV one-loop 1PI amplitude.
2. The complete native UV counterterms and missing hard one-loop topology
   inventory have **not** been independently derived.
3. The signed `GammaCC`/Majorana-Weyl normalisation is not certified.
4. Missing native scalar operators / full physical scalar-threshold matching
   remain separately uncertified.
5. No claim of a complete physical `MSbar` Wilson coefficient or a physical
   neutrino-mass prediction is made. Native pole residuals must **not** be set
   to zero numerically or absorbed by a fitted coefficient.

The resulting audit is a concrete route to identifying which native PL
source sectors require additional UV/EFT1 diagram calculations. It is not a
replacement for those calculations.

## Files in this ZIP (relative to project root)

| Path | Action | Purpose |
|---|---|---|
| `Numerical/parallel_t3/wilson/native_pole_closure_batch7.py` | Add | Non-destructive symmetry decomposition, native-source grouping, CT/B0 checks, pole evidence and certification gates |
| `scripts/RunT3PoleClosureBatch7.py` | Add | All-16 source-resolved runner, saved-source hashes, scale scans, conditional SMEFT evolution and output ledgers |
| `scripts/RunT3PoleClosureBatch7.ps1` | Add | PowerShell launcher |
| `tests/test_t3_pole_closure_batch7.py` | Add | Isolated regression tests, including mixed-source `B_-1` failure mode |
| `README_T3_NATIVE_POLE_CLOSURE_BATCH7.md` | Add | This document |

## Install / run (PowerShell)

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "$HOME\Downloads\t3_native_pole_closure_batch7_all16.zip" `
  -DestinationPath . -Force

.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_pole_closure_batch7.py

# First: diagnose and recover B_-1 only.
.\.venv\Scripts\python.exe scripts\RunT3PoleClosureBatch7.py --model B_-1

# Then: all 16 models, including conditional low-energy running.
.\.venv\Scripts\python.exe scripts\RunT3PoleClosureBatch7.py --resume

# Inspect all16 summary.
.\.venv\Scripts\python.exe scripts\RunT3PoleClosureBatch7.py --show
```

Or use the wrapper:

```powershell
.\scripts\RunT3PoleClosureBatch7.ps1 -Model B_-1
.\scripts\RunT3PoleClosureBatch7.ps1 -Resume
.\scripts\RunT3PoleClosureBatch7.ps1 -Show
```

The B_-1 first command may reveal a further physical consistency blocker;
there is no fabricated success fallback. `--resume` reruns failed/missing models.
To regenerate an existing result with a changed low-energy scale or from a
pole-only `--skip-low-energy-running` run, use `--resume` with default running
or `--overwrite` (which **only** removes the Batch-7 directory).

## Output contract

`output/parallel_t3/native_pole_closure_all16_batch7_v1/`

- `summary.json`: all16 completion/blockers, source metrics, individual paths.
- `all16_source_pole_closure.csv`: compact 16-model comparison.
- `<model-tag>/source_pole_and_matching_ledger.json`: detailed raw source
  symmetry diagnostics, native-versus-EFT1 residuals, pole-normalisation
  *diagnostic* fits and explicit matching-completeness gates.
- `<model-tag>/source_poles_and_counterterms.npz`: all raw complex source
  arrays, their totals and difference matrices.
- `<model-tag>/matching_scale_scenario_matrices.npz`: branch Wilson tensors.
- `<model-tag>/conditional_mass_and_mixing_scan.csv`: 9 conditional branch
  results when low-energy running is enabled.

Existing Batch-6 `15/16` outputs remain intact and are cross-checked when
present. If `B_-1` finishes in Batch 7 it is a **completed conditional numerical
audit**, not a newly established physical one-loop matching coefficient.

## Validation status

Isolated pure-NumPy tests can be executed without Windows source outputs.
The true all16 numerical and repository integration run is available only on
the local project with the saved Batches 1–5 data. Do **not** interpret
isolated unit tests as evidence of successful all16 execution.

## Batch 7 resume hotfix (v2)

The original Batch 7 runner raised `NameError: name 'skip_running' is not defined`
when `--resume` evaluated an existing completed record. The full replacement
`scripts/RunT3PoleClosureBatch7.py` now reads the parsed command-line flag
`args.skip_low_energy_running`. The saved `B_-1` result is reusable, and the
remaining models run without deleting Batch 7 outputs. The fix does not change
the physics calculations or any matching inputs.

From the project root, after installing this ZIP:

```powershell
.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_pole_closure_batch7.py tests\test_t3_batch7_resume_cli.py
.\.venv\Scripts\python.exe scripts\RunT3PoleClosureBatch7.py --resume
.\.venv\Scripts\python.exe scripts\RunT3PoleClosureBatch7.py --show
```

This workstream's practical completion condition is a reproducible full run
with 16/16 recorded outputs, a saved input/output ledger and no traceback.
The native-pole discrepancy remains a diagnostic in the output, not a reason
to withhold working numerical scenario results.
