# Section 5.3 — Zee-style compensated sensitivity with individually fitted general T3 inputs

## Installation and changes

This is a **full replacement** of three existing project files (not a patch):

- `scripts/RunT3GeneralZeeSensitivity.py` — numerical driver and merged-fit source checks;
- `tests/test_t3_general_zee_sensitivity.py` — pure source-selection/transformation contracts;
- `docs/T3_GENERAL_ZEE_SENSITIVITY.md` — this guide.

Extract the archive into the **local repository root**, preserving the directory structure. No change or push to GitHub is performed. The old fitted inputs and existing scan outputs are read-only.

## Verified benchmark selection (GitHub main, 4 November 2026)

The default `--fit-run-id merged` loads the saved accepted-candidate manifest:

`output/studies/general_scalar_extended/thesis_general_all16/continue_general16_v1/merged_all16.json`

The manifest reports **16 accepted and zero blocked** under the exploratory `chi2 <= 5.86` fit screen. It is assembled from:

- `OriginalGeneralCampaign`: 11 model-specific fits from `thesis16_general_c5_full_v1`;
- `GeneralFitContinuation`: improved A(-2), B(+1), C(-1), D(0), D(+2) fits from `continue_general16_v1`.

For example, A(-2) must use the continuation (`chi2 = 5.520689216404163e-11`), **not** the original (`chi2 = 34.52948208294043`). C(-1) uses its improved continuation (`chi2 = 1.79419851843013e-12`). The very small `2.7451024734292607e-12` E(+2) result comes from the original campaign and is already included in the merged manifest.

Each selected model's `result_path` must equal its expected location within these **two allowed directories**. Its `fitted_config.json` and any `general_uv_overrides.json` are always read from that same run and model folder. The driver refuses malformed/incomplete manifests, poor fits, model mismatches, mismatched saved chi-square, inconsistent saved 18-entry `y1/y2` optimiser parameters, inconsistent general UV overrides, and unexpected EFT-matching certification. It saves the per-model source selection and source SHA-256 hashes with results.

**Important:** The existing fit optimised 18 real entries of `y1` and `y2`. It did **not** independently optimise every scalar coupling or each additional conjugate Yukawa. The latter retain the saved general-boundary assumptions.

### Source compatibility

The continuation `result.json` format differs from the original campaign: the nonmixed continuation stores its general-RGE/matching provenance in `matching_and_rge_provenance`, and C(-1)'s independent quartic/bilinear fractions are nested in that record rather than the original top-level keys. The driver accepts each verified format and does not mix their parameter values.

`--fit-run-id thesis16_general_c5_full_v1` is retained **only** as an explicit option to reproduce the old campaign. If omitted, the driver uses `merged`.

## Physics calculation

The three UV transformations remain:

- **T1:** `lambdaT3 -> gamma*lambdaT3`, `y1,y2 -> gamma**(-1/2)*(y1,y2)`;
- **T2:** `lambdaT3 -> gamma*lambdaT3`, `y2 -> gamma**(-1)*y2`;
- **T3:** `y1 -> gamma*y1`, `y2 -> gamma**(-1)*y2`, `lambdaT3` unchanged.

The nominal tensor `lambdaT3*(y1@y2.T + y2@y1.T)` is invariant at the UV scale. Other general-sector couplings, including independent extra quartics, allowed scalar bilinears, and conjugate Yukawas, are held fixed at the UV boundary. For B/C, the numerical adapter's independent quartic input fraction is compensated to prevent it from changing when `lambdaT3` changes.

Each gamma point re-evaluates the existing general one-loop UV RGE from `1e7` to `1e5` GeV, performs the existing fermion threshold projection, evolves the general scalar EFT between `1e5` and `1e3` GeV, evaluates the available model-specific general **one-scale candidate C5**, then evolves the SM+Weinberg theory to `1e2` GeV.

This is **not certified complete sequential general EFT matching**. The original general Matchete C5 is conditional for 14 nonmixed models, and B(-1)/C(-1) use separate conditional mixed-mass one-scale matching with an uncertified signed quartic map. Full fermion-first Wilson matching and scalar mass-basis threshold matching are not implemented here. The scan never replaces missing general matching with historical restricted kernels.

## Commands (PowerShell)

From the **local** `T3neutrinomassmodelpipeline` project root:

```powershell
# If the archive is in Downloads, change the path to its actual location.
Expand-Archive -LiteralPath "$HOME\Downloads\T3_General_Zee_MergedFits_Update.zip" `
    -DestinationPath . -Force

# Pure unit tests: no RGE integrations or Wolfram calls.
python -X utf8 -m pytest -q tests\test_t3_general_zee_sensitivity.py

# Verify exactly which general fit A(-2), B(-1), C(-1) will use.
python -X utf8 -m scripts.RunT3GeneralZeeSensitivity `
    --models A_-2,B_-1,C_-1 --preflight

# First, repeat A(-2) with IMPROVED continuation Yukawas.
# Use a new run ID: the previous zee_general_a2_v1 contains the OLD fit.
python -u -X utf8 -m scripts.RunT3GeneralZeeSensitivity `
    --models A_-2 --run-id zee_general_a2_merged_v1 `
    --gamma-min 0.99 --gamma-max 1.01 --points 11

# If gamma=1 reproduces the continued saved fit, scan all 16.
python -u -X utf8 -m scripts.RunT3GeneralZeeSensitivity `
    --models all --run-id zee_general_all16_merged_v1 `
    --gamma-min 0.99 --gamma-max 1.01 --points 11
```

The `--preflight` command checks the merged selection and actual general sources, **without numerical integration or output creation**. The full run aborts rather than overwriting an existing `--run-id`; give each attempt a unique ID.

## Outputs

Under

`output/studies/general_scalar_extended/zee_style_general_sensitivity/<run-id>/`

- `summary.json`: model status, selected source and saved chi-square, manifest SHA-256;
- `points.csv`: direction, gamma, five oscillation observables and chi-square;
- `<model_tag>/scan.json`: every evaluated point, the exact source run (original or continued), SHA-256 source hashes, matching limitations and gamma=1 verification;
- `<model_tag>/solar_sensitivity.png`: the **T1/T2/T3 trajectories only**, if matplotlib is available.

A scan requires **all three gamma=1 predictions to reproduce the selected model's saved values**. If this fails, the model is recorded as `Blocked`; that is a diagnostic, not permission to substitute an older fit. Failed noncentral points are saved with errors, not interpolated.

### Subsequent NuFIT/JUNO sensitivity ovals

The numerical driver currently generates trajectories but **does not yet overlay** NuFIT or JUNO confidence regions. Build that separate plotting layer from the *new merged-fit* `points.csv`, not historical reduced-potential sensitivity files. Existing repository target: `configs/nufit_6_1_ic24_no.json` is a diagonal Gaussianised **NuFIT 6.1** approximation, not the official correlated likelihood. Existing `docs/FROZEN_LAMBDA_TARGET_REFERENCE.json` includes a **JUNO six-year projection** for solar parameters; it is a future precision estimate, not a measurement. Distinguish (a) model sample-covariance sensitivity ovals, (b) NuFIT approximate confidence regions, and (c) JUNO projected precision. No official experimental 2D contour is implied by independent 1-sigma inputs.

## Validation scope

Pure Python tests verify 21 contracts, including all three scaling directions, nominal tensor invariance, independent scalar boundary invariance, mixed B/C compensation, merged source provenance/path selection, prevention of poor-fit fallback, continuation B/C scalar input reading, and rejection of mismatched fitted Yukawas. They do **not** execute the full RGEs or Wolfram/Matchete in this artifact environment. A real local `--preflight` and gamma=1 scan remain necessary.
