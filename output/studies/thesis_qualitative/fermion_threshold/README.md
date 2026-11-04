# General T3 physical fermion-threshold study

Code: `scripts/RunT3GeneralFermionThresholdScan.py`.

**November 2026 correction (full-file replacement):** passes all four required keyword-only tolerances to the production `_heavy_mass_basis_data` helper and reads mixed-scalar candidate fractions from the actual saved fit format. Original B(-1) records use top-level `extra_quartic_fraction` / `extra_mixed_scalar_bilinear_fraction`; the C(-1) continuation uses `matching_and_rge_provenance.quartic_fraction` / `.mu2_fraction`. The driver now checks this provenance during preflight and rejects missing/conflicting or modified scalar inputs. These are driver corrections; no physics kernels have changed. Previous outputs are retained unchanged.

Commit inspected: `6ebe9937a72880c6f180bcd026261779800ac71f` (`main`).

## Study definition

- Sixteen **manifest-selected** general-Lagrangian fitted configurations from `output/studies/general_scalar_extended/thesis_general_all16/continue_general16_v1/merged_all16.json`.
- Three physical UV heavy-fermion mass matrices: `MF(UV) = r MF_0(UV)` for `r = 0.5, 1, 2`. The *fermion decoupling scale* changes simultaneously: `muF(r) = r muF_0`. The UV Yukawa matrices, scalar masses, quartics and UV scale are not refitted or changed.
- Strict hierarchy: UV with all heavy fields -> integrate out **all** fermions at the varied `muF` -> scalar-only EFT1 containing **both** heavy scalars -> integrate out **all** scalars at fixed `muS = 1000 GeV` -> SM + Weinberg operator evolved to `100 GeV`.
- Baseline at `r=1` is evaluated **first**, and compared with the five saved fit observables. By default each relative difference must be <= 0.5%; a failed baseline prevents the model's varied-mass points from being calculated.
- A/D/E and other non-mixed models: actual fitted `general_uv_overrides.json`, `GeneralRunningBenchmark.initial_state`, `UVRunner`, `IntermediateScalarRunner`, and **full-general** `GeneralT3FlavourC5`. The generated coefficient and canonical export fingerprints are checked. No restricted C5 coefficient is substituted.
- B(-1)/C(-1): fit-selected `fitted_config.json`, native quadratic+quartic candidate RGE payload, projected scalar-only state, evolved holomorphic scalar mixing, mixed-scalar mass eigenvectors, existing conditional one-scale Matchete coefficient. The fixed-`1e5`-GeV evaluator is *not* called; two native integration intervals explicitly use the varied `muF`.
- Every completed result saves the UV physical `MF`, evolved `MF(muF)`, heavy-fermion masses, both scalar masses (mixed eigenmasses for B/C), threshold endpoints, full C5 at scalar threshold and after SMEFT running, neutrino mass eigenvalues, five oscillation observables, and relative/directional changes.

## Scope limitation (essential)

The matched C5 is a **conditional one-scale** general Matchete coefficient at the scalar threshold. The computation respects the fermion-first *running and field-content order* but is **not** a fully certified EFT1 Wilson-coefficient matching calculation. In B/C the two direct quartics use the existing **uncertified** sign convention `lambda_Matchete = -lambda_RGBeta` and native mixed-scalar candidate basis. This is suitable for a *conditional qualitative study*, not precision physical predictions. Model-dependent failures are `Blocked` and are never replaced by restricted kernels.

## PowerShell from repository root

```powershell
python -X utf8 -m pytest -q tests/python/test_general_fermion_threshold_scan.py
python -X utf8 -m scripts.RunT3GeneralFermionThresholdScan --preflight --run-id preflight_fixed_v2 --models A_-4,B_-1
python -X utf8 -m scripts.RunT3GeneralFermionThresholdScan --run-id pilot_fixed_v2 --models A_-4,B_-1
python -X utf8 -m scripts.RunT3GeneralFermionThresholdScan --preflight --run-id preflight_all16_fixed_v2 --models all
python -X utf8 -m scripts.RunT3GeneralFermionThresholdScan --run-id all16_fixed_v2 --models all
```

Use a **new** `--run-id` for each invocation. The driver refuses to overwrite previous study runs.

The first two pilot entries may require installed Python physics dependencies and the `wolframscript` executable (non-mixed full-general channel matching). Native B/C conditional one-scale evaluation is performed by existing Python/mPMath code.

Each run produces:

- `<run-id>/settings.json` with SHA and scale/factor metadata.
- `<run-id>/<model-tag>/preflight.json`, `factor_1p0.json`, `factor_0p5.json`, `factor_2p0.json`, or explicit `blocked.json`.
- `<run-id>/summary.json`, `summary.csv`, `trends.csv`, `representation_dependence.json`.

`representation_dependence.json` aggregates **only actually computed models**. Family medians combine distinct fitted UV Yukawa textures and are **not** isolated Clebsch-Gordan / representation coefficients.

## Validation in the file-authoring environment

`16 passed` for isolated study contract tests. These cover physical mass rescaling, unchanged fitted Yukawas, scale hierarchy, baseline gating, path safety, ratio/direction calculations and regression guards for the fixed-threshold evaluator. The source repository itself is **not present in the execution container**, so full RGBeta/Matchete runs of the **corrected** driver have not been executed here. The user reported a successful B(-1) three-point scan with the previous driver (m3 = 0.08740539, 0.05246321, 0.03054945 eV at r = 0.5, 1, 2 respectively); these numbers were not independently rerun in this workspace. The A(-4) and remaining nonmixed models were blocked by the missing keywords, and C(-1) by the fit-format mismatch. The corrected run must be completed locally before claiming all-16 trends.
