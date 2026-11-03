# Thesis 5.1 — shared A(-2) fitted general-Lagrangian Yukawas on all 16 T3 models

## Objective

Reuse the **saved full-general A(-2) fitted** 3x3 real non-diagonal `y1`, `y2`
*unchanged at the UV boundary for every one of the 16 ordinary T3 models.*
Run each target model's canonical **full-general UV renormalisable RGEs** from
`1e7 -> 1e5 GeV`; integrate out the **fermions**, then its **general scalar-only
EFT1 RGEs** from `1e5 -> 1e3 GeV`; evaluate existing conditional one-scale
`C5`; run SM+Weinberg `1e3 -> 1e2 GeV`; compute neutrino masses and the five
oscillation observables. No new fitting, no historical results overwritten,
no GitHub changes.

Sources (read-only):

- `output/studies/general_scalar_extended/thesis_general_all16/thesis16_general_c5_full_v1/A_m2/fitted_config.json`
- `output/studies/general_scalar_extended/thesis_general_all16/thesis16_general_c5_full_v1/A_m2/general_uv_overrides.json` (fixed A(-2) `y2Conj`)

The source A(-2) result has `chi2_final=34.52948208294043`, **a poor fit** under
the campaign's `chi2 <= 5.86` screening convention. This run is a controlled
**fixed-input comparison**, not a survey of viable independently fitted models.

`y1`, `y2` are copied **exactly** from the general A(-2) fitted config.
Additional `y1Conj` / `y2Conj` Yukawas exist only for representations allowing
them. The exact **saved A(-2) `y2Conj` matrix** is copied into any target
model where `y2Conj` is allowed (A(-2), D(-2), E(-2)). In the original A(-2)
general run, `y2Conj` was a fixed auxiliary Yukawa equal to 0.05 times its
*pre-fit* reference `y2`; it is **not** 0.05 times the newly fitted A(-2) `y2`.
For models permitting `y1Conj` instead (A(0), D(0), E(0)), that interaction
does not exist in A(-2). It is assigned the declared candidate convention
`y1Conj = 0.05 * shared A(-2) y1`. No extra Yukawa is independently refitted.
B/C models have no additional conjugate-scalar Yukawa.

All other initial conditions remain model-specific, loaded through the same
`Numerical/benchmarks/fitted_reference/` reference configurations and new
candidate scalar interactions used by the general-model fitter, with scalar
additional quartics `0.1*abs(lambdaT3)` and mass bilinears
`0.01*min(mS1Sq,mS2Sq)` by default. These fixed model-specific values do not
represent a fixed *full* UV input point across representations; only `y1,y2`
are identical. Class E retains its reference topology-coupling convention.

## Installed file

- `scripts/RunT3GeneralSharedA2.py`: new, complete standalone script;
  reuses current repository numerical kernels. It does not modify production.

## Commands (PowerShell, from the local project root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

# Install the package at the repository root.
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\T3_general_shared_A2_all16_runner.zip" `
  -DestinationPath . -Force

# Check local Python dependencies and source compatibility.
python -X utf8 -m scripts.RunT3GeneralSharedA2 --self-test

# Verify each model's general betas, Matchete C5, B/C mixed-scalar inputs.
python -X utf8 -m scripts.RunT3GeneralSharedA2 `
  --run-id thesis51_shared_a2_v1 --preflight

# Calculate all 16 models once; no Sobol optimisation.
python -u -X utf8 -m scripts.RunT3GeneralSharedA2 `
  --run-id thesis51_shared_a2_v1

# If some models fail and their local prerequisites have been repaired:
python -u -X utf8 -m scripts.RunT3GeneralSharedA2 `
  --run-id thesis51_shared_a2_v1 --resume

Import-Csv .\output\studies\general_scalar_extended\shared_general_A_m2\thesis51_shared_a2_v1\summary.csv | `
  Format-Table model,status,chi2,k_to_A_m2_real,epsilon_to_A_m2,m1_ev,m2_ev,m3_ev -AutoSize
```

The runner uses `--canonical-run-id full16_general_v1`. It needs the already
saved canonical general beta exports, compiled/source Matchete C5 expressions,
Wolfram-dependent *previously exported* files, and the B/C candidate beta and
one-scale mixed-scalar matching sources. `--preflight` does not integrate or
optimise; it reports per-model missing prerequisites. The default run does
not invoke Wolfram or regenerate missing data. Missing or unsupported models
are recorded as `Blocked`, **not** replaced by reduced-potential C5.

## Outputs (new-only)

`output/studies/general_scalar_extended/shared_general_A_m2/thesis51_shared_a2_v1/`

- `settings.json` — frozen input and source hashes
- `summary.csv` / `summary.json` — per-model numerical status, C5 norms,
  mass spectra, mass splittings, mixing observables, chi2
- `<model_tag>/status.json` — per-model success/failure and provenance
- `<model_tag>/result.json` — full C5 matrices at scalar and low scales
- `fig_shared_a2_general_chi2.png/.pdf` — conditional fit-quality comparison
- `fig_shared_a2_general_C5_residual.png/.pdf` — complex C5 matrix-shape
  residual relative to A(-2), if A(-2) matched successfully

With A(-2) as reference,

`k = Tr(C5_A2^dag C5_i) / Tr(C5_A2^dag C5_A2)`,  
`epsilon = ||C5_i - k*C5_A2||_F / ||C5_i||_F`.

Masses are extracted as Takagi singular values using
`m_nu = -(246.22 GeV)^2 * C5/2`, converted GeV to eV. All physically relevant
mixing observables are obtained via the existing `low_energy_fit` routine.

## Scientific scope

The 14 nonmixed models use the existing **canonical full-general UV/EFT1
runner and general Matchete one-scale C5**. The B(-1)/C(-1) mixed-scalar cases
use a **different existing conditional** candidate scalar-mass-basis and
one-scale Matchete matching adapter because the generic full-general pipeline
is not certified for them. This is an explicitly method-dependent exploratory
comparison; the cross-method normalisation is not certified. Results may be
blocked by unavailable candidate beta files or missing basis transformations.

**Not performed or certified:** full physical Matchete/RGBeta parameter map,
complete fermion EFT Wilson-operator matching/running, or full physical
scalar mass-eigenstate threshold matching for every model. Therefore do not
label the 16 outputs 'fully matched UV -> SMEFT predictions', even if all
conditional evaluations succeed.

**No independent model refits** are performed: a poor fixed-input chi2 means
the same A(-2) Yukawas fail to describe that representation at this benchmark;
it does not exclude that model after its Yukawas are refitted.
