# Batch 7 — B/C zero-additional-coupling numerical benchmark

## Scope

This is a **new independent numerical run** in the **existing, uncertified native RGBeta scalar basis**, using the repository's saved UV and EFT1 one-loop beta exports. It does **not** replace the production Lagrangian, original fitted benchmark, Stage30–36 outputs, or physical threshold matching.

The current public `main` was inspected at commit `a02e04c6def5cf8522b2f612ccd6bff7d9963b91`. Supported cases: `B_-1`, `C_-1` only. Other models were covered in Batches 1–6; they have no direct additional quartic contribution to the selected one-loop C5 graph.

## Inputs and boundary condition

At `mu_UV = 10^7 GeV`, the original fitted `B_-1` or `C_-1` **ordinary** numerical state is loaded without refitting. Its two scalar diagonal masses, SM parameters, fermion masses, Yukawas and original `lambdaT3` remain the saved fitted benchmark values. Only the model-specific additional scalar parameters are set to zero in memory:

- Six additional quartics already present in the existing RGBeta state: `lambdaHHdagS2S2`, `lambdaHHdagS1barS1bar`, `lambdaHHdagS1barS2barCross`, `lambdaS1bar2S2bar2`, `lambdaS1barS2S2bar2`, `lambdaS1S1bar2S2bar`.
- One new candidate quartic: `lambdaT3GeneralQ110101I1`.
- The allowed mixed quadratic: `muT3GeneralS1S2` (dimension two).

Thus **all seven additional quartic model-slots and the one extra quadratic are exactly zero**, not just numerically small. Their ordinary scalar mass parameters and original potential quartics are **not** set to zero.

### Two-stage evolution

UV: `10^7 -> 10^5 GeV` with the existing candidate **UV** RGBeta export containing the full ordinary + seven additional quartic betas + the extra quadratic beta. A separate saved Stage31 UV export is used to check that its original ordinary/general equations have equal numerical derivatives at the zero-`mu2` boundary.

At `mu_F = 10^5 GeV`, the existing `project_general_after_fermion` tree-level scalar boundary is used. All additional quartic and quadratic values are explicitly transported without modification; **no finite fermion-threshold Wilson or scalar correction is calculated**.

EFT1: `10^5 -> 10^3 GeV` with the existing candidate **scalar-only EFT1** RGBeta export. Both intervals use `solve_ivp(method='DOP853', rtol=1e-8, atol=1e-11, max_step=0.35 in ln(mu))` and source-key inventories are checked before integration. Source-file SHA-256 fingerprints are written alongside the outputs.

At `mu_S = 10^3 GeV`, record the generated quartics in RGBeta coordinates and `|mu2|/min(mS1²,mS2²)`. The ratios `|lambda_direct|/|lambdaT3|` are coupling-coordinate diagnostics, **not C5 ratios**.

## Perturbative interpretation

With UV additional quartics zero, saved one-loop RGBeta sources regenerate three B/C additional quartics. Two are the operators that are direct candidates in the saved Matchete one-loop C5:

- `lambdaHHdagS2S2` (candidate correspondence `lambdaT3ExtraQ020020I1`),
- `lambdaHHdagS1barS1bar` (candidate correspondence `lambdaT3ExtraQ020200I1`).

The *numeric equality* of RGBeta and Matchete coupling values is **not verified**. At leading logarithmic order their generation is schematically

`lambda_extra(mu_F) ~ [beta_extra^(1)(mu_UV)/(16*pi²)] log(mu_F/mu_UV)`.

For `100:1` scale ratios, `ln(100)/(16*pi²) ~ 0.02916`. The actual source coefficients, signs, cancellation and running matter; the script integrates the full coupled candidate one-loop equations instead of multiplying by a guessed coefficient.

If a generated additional quartic is inserted into a genuinely one-loop C5 amplitude, its contribution begins formally one RG-loop order beyond the fixed-order one-loop baseline. This makes a strict fixed-order leading-one-loop truncation defensible under the declared UV boundary; it is **not** a bound on an RG-resummed or numerically fitted observable.

## Optional `--c5-preview`

The optional preview calls the **already existing** `RunGeneralT3FrozenScalarC5.build_scenario` and `RunGeneralT3NumericalC5.numerical_model` routines, loading the saved direct one-scale Matchete coefficient. It evaluates **two one-scale scenarios with exactly the same inputs except the two direct quartics**:

1. Baseline: `lambdaT3` from the current EFT1 candidate endpoint, two direct additional Matchete quartics zero.
2. Preview: same `lambdaT3`, with `lambdaT3ExtraQ020020I1 = lambdaHHdagS2S2` and `lambdaT3ExtraQ020200I1 = lambdaHHdagS1barS1bar` **numerically identified by assumption**.

Both scenarios use frozen original Yukawa matrices, fermion masses and **frozen reference scalar diagonal masses**, not running mass eigenstates. The preview refuses to evaluate if the computed quadratic mixing is appreciable, because the frozen-diagonal one-scale matching expression cannot represent nonzero general scalar mixing.

The reported sensitivity is

`epsilon_C5 = ||C5_identity_map - C5_zero_extra||_F / ||C5_zero_extra||_F`.

This is an **uncertified sensitivity preview**: no Matchete/RGBeta tensor normalization certificate, no running Yukawa or scalar eigenstate inputs in the C5 evaluation, no Wilson running or complete two-threshold MSbar matching. **Do not present epsilon_C5 as a measured physical neutrino-mass error bar**, even if small.

## Installation and commands (PowerShell)

Copy all files from the ZIP into the repository root, preserving `scripts/`, `tests/python/`, and `docs/`:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Scalar_Relevance_Batch7.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_zero_extra_boundary_batch7.py
python -m scripts.RunT3ZeroExtraBoundaryBenchmark --all
Get-Content -Encoding UTF8 output\scalar_relevance\zero_boundary_benchmark\summary.md
```

**Optional exploratory C5 ratio** (requires the existing frozen-reference Matchete one-scale coefficient and mass-basis probe files):

```powershell
python -m scripts.RunT3ZeroExtraBoundaryBenchmark --all --c5-preview
Get-Content -Encoding UTF8 output\scalar_relevance\zero_boundary_benchmark\summary.md
```

Outputs, all generated **locally**:

```text
output/scalar_relevance/zero_boundary_benchmark/
    summary.json
    summary.md
    B_m1/result.json
    C_m1/result.json
    B_m1/one_scale_preview/baseline/B_m1_final_numerical_c5.json        # optional
    B_m1/one_scale_preview/identity_map/B_m1_final_numerical_c5.json    # optional
    C_m1/one_scale_preview/baseline/C_m1_final_numerical_c5.json        # optional
    C_m1/one_scale_preview/identity_map/C_m1_final_numerical_c5.json    # optional
```

The numerical C5 helper also writes one `*_mnu_eV.csv` beside each optional C5 JSON output.

## What is and is not verified

- **Verified from the repository:** supported B/C operator inventory, both beta-export schemas, the existing Stage31 numerical UV state builder, the existing scalar-only field projection, and the optional frozen one-scale Matchete numerical C5 API.
- **Implemented in this batch:** independent zero-extra UV and EFT1 candidate coupled integrations, zero-boundary and stage-continuity checks, independent initial UV beta-system cross-check, provenance hashes, per-model generated-coupling magnitudes, optional identity-map C5 sensitivity.
- **Tests run here:** offline schema guards, integration kernel, report and C5 norm tests only. Actual B/C integration and optional Matchete preview **must run on the user's checkout**.
- **Not verified:** physical RGBeta–Matchete component/intertwiner normalization, true UV-to-EFT1 finite matching, full sequential Weinberg coefficient, generalized mass mixing at one loop, physical neutrino mass or phenomenological acceptance.

The batch is intended to **stop** after producing the B/C zero-extra-coupling benchmark and optional sensitivity. No thesis-level scope expansion is implied.
