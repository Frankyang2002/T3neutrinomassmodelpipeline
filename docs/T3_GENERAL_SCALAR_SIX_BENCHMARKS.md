# Six general-scalar benchmarks: reuse existing 16 × 3 Sobol campaign

## Scope

**Already implemented in the repository:** `scripts/RunT3GeneralScalar48Campaign.py` is the historical general-quartic exploratory *Sobol + local least-squares* runner for all 16 ordinary T3 models and three scalar-topology/UV scenarios. Do not start over with `Numerical/fitting/SobolBenchmarkSearch._build_evaluator` because that entry point uses the earlier reduced-potential `FinalC5TrajectoryEvaluator`.

The three actual historical scenario IDs (from `scripts/FitT3BPlus1GeneralThreeStudies.py`) are:

| Scenario | Topology lambdaT3 at UV | UV RG |
|---|---:|---|
| `large_uv_frozen_0p5` | 0.5 (E: 0.25) | Frozen; boundary imposed at 100 TeV |
| `small_uv_frozen_5e-9` | 5e-9 (E: 2.5e-9) | Frozen; boundary imposed at 100 TeV |
| `small_uv_running_5e-9` | 5e-9 starting value (E: 2.5e-9) | Coupled 1e7 GeV to 1e5 GeV |

The two Yukawa strategies are **individually optimised** and **fixed A(-2) Yukawas**. Each scenario has 16 models × 2 strategies, so six configurations and 96 comparison cells. Only the 48 individually optimised cells are part of the original bulk campaign; the new local runner attempts the 48 fixed-A(-2) comparisons with the corresponding existing general-scalar candidate evaluators. No historical results are changed.

**IMPORTANT numerical result:** The currently inspected saved `bulk48_v1/campaign.json` contains 48 cells and 20 conditional acceptances. The A(-2) saved donor is accepted in `large_uv_frozen_0p5` (chi2 ≈ 1.3288e-4), poor in `small_uv_frozen_5e-9` (chi2 ≈ 2.1368e12), and has no mass-safe C5 evaluation for `small_uv_running_5e-9`. The new runner **does not invent a missing A(-2) donor**. Thus **96 valid numerical evaluations cannot yet be claimed**.

## Files added

- `scripts/CompareT3GeneralScalarSharedA2.py` — read-only loader for existing 48 Sobol results; optional calculation of shared-Yukawa comparisons using the existing model-specific general-scalar candidate evaluators.
- `tests/test_t3_general_scalar_shared_a2.py` — synthetic-source tests of 48/96 coverage, source status checks and missing donors.
- `docs/T3_GENERAL_SCALAR_SIX_BENCHMARKS.md` — this guide.

No existing scripts, `continuous_v1` data, original lambda display code, or historical results are replaced. No GitHub write.

## PowerShell commands

Install the ZIP into `C:\thesistemp\T3neutrinomassmodelpipeline` preserving folders, then:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m pytest -q tests\test_t3_general_scalar_shared_a2.py

# Check 48 saved Sobol fit cells and inspect A(-2) donor feasibility.
python -X utf8 -m scripts.CompareT3GeneralScalarSharedA2 --campaign bulk48_v1 --check

# Construct 96-row table WITHOUT any new physics evaluations.
python -X utf8 -m scripts.CompareT3GeneralScalarSharedA2 `
    --campaign bulk48_v1 --run-id shared_a2_dry_v1 --dry-run

# Launch genuinely new shared-Yukawa calculations, only where an accepted
# scenario-specific A(-2) fit exists. This is the safe default.
python -u -X utf8 -m scripts.CompareT3GeneralScalarSharedA2 `
    --campaign bulk48_v1 --run-id shared_a2_v1
```

The optional command below additionally attempts the frozen-low-scalar scenario using the **poor-fit A(-2) donor**, with that limitation explicitly recorded in every output row:

```powershell
python -u -X utf8 -m scripts.CompareT3GeneralScalarSharedA2 `
    --campaign bulk48_v1 --run-id shared_a2_poor_donor_v1 --allow-poor-donor
```

Choose a new run ID when changing settings. Use `--resume` **only with identical settings** after interruption. Avoid re-running historical fits until their status/provenance has been reviewed.

## Outputs

All new files live under:

`output/studies/general_scalar_extended/six_benchmarks/<run-id>/`

- `a_m2_donors.json`: same-scenario A(-2) Yukawa source files, hashes, chi2 and donor eligibility.
- `shared_results.json`: saved/resumable per-model shared candidate calculation status and observables, or honest unavailable/failed labels.
- `comparison.json`: exactly **96 rows** (48 previously saved optimisation outcomes, 48 shared cases), retaining unavailable rows explicitly and documenting physics limitations.
- `comparison.csv`: sortable external CSV with model, scenario, mode, chi2, available oscillation observables and source status.

This is not a production full-scalar physical calculation. The 48-campaign's A/D/E extra holomorphic quadratic interactions were omitted. B(-1)/C(-1) use conditional one-scale mixed-scalar Matchete C5; B(+1) uses its Batch21 source; other models use the Bulk48 hierarchical/legacy source. None of these is a certified full general-potential sequential EFT calculation. B/C shared comparisons use the original candidate `quartic_fraction=0.1, mu2_fraction=0` by default and require separate screening; they are not assumed to have the identical scalar boundary as the 14 bulk models.

The original Display Lambda RGE and T1/T2/T3 sensitivity plotting code is not edited. **This package does not generate genuine general-scalar gamma sensitivity scans or ovals.** Those require a full evaluator per fitted point and T1/T2/T3 rescans; do not copy historical reduced-potential curves.

## Next numerical decision

Run `--check` and inspect donor statuses. The `small_uv_running_5e-9` A(-2) donor requires a new mass-safe Sobol/fit search; simply calling a shared comparison without a donor would not constitute the sixth benchmark. The current source flags full sequential general-scalar matching as uncertified, so any further fits must retain their conditional-approximation labels.
