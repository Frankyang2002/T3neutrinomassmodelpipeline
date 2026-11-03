# Section 5.1: dual fixed-Yukawa general-T3 representation comparison

## Purpose

Run **two separate, fixed-input scans** across the 16 ordinary T3 models, using the previously fitted general-Lagrangian **D(-2)** and **A(-4)** Yukawa textures. No new fitting, no modifications to the repository's production physics code, no GitHub changes, and no restricted-Lagrangian substitutions.

Sources (read-only, canonical general-campaign fits):

- `output/studies/general_scalar_extended/thesis_general_all16/thesis16_general_c5_full_v1/D_m2/fitted_config.json`
- `output/studies/general_scalar_extended/thesis_general_all16/thesis16_general_c5_full_v1/D_m2/general_uv_overrides.json`
- `output/studies/general_scalar_extended/thesis_general_all16/thesis16_general_c5_full_v1/D_m2/result.json`
- The corresponding `A_m4/` files.

Original candidate fit chi2: D(-2) = **4.156852052877585e-14**, A(-4) = **2.7184618585331287e-7**. These are conditional one-scale general-matching results, not fully certified sequential SMEFT predictions.

For each source model `s`, fix the **18 real numbers** in `Y1_s,Y2_s` identically at `mu_uv=1e7 GeV` for every target representation `i`. The target model's other reference scalar/SM/heavy masses and general interactions are preserved, using the same conventions as `scripts.RunT3GeneralAll16FromReference`. In particular, the allowed extra quartics are initially `0.1|lambdaT3|`, and the allowed extra bilinear is initially `0.01 min(mS1Sq,mS2Sq)`. Target models therefore **do not have identical complete UV parameter points**; only their transferred `Y1,Y2` are identical within each sweep.

**Extra Yukawa convention:** For the D(-2) sweep, the actual saved D(-2) additional `y2Conj` matrix is used wherever that Yukawa is allowed, preserving the D(-2) source fit. If a target admits `y1Conj` instead, set it to `0.05 Y1_D` (not a direct transfer of a source parameter). For the A(-4) sweep (which has no additional Yukawa), any target `y1Conj`/`y2Conj` is assigned `0.05 Y1_A`/`0.05 Y2_A` respectively. This introduces a **different auxiliary-coupling prescription** between the two sweeps; interpret small texture-dependent differences with this qualification.

## Numerical pipeline and limitation

- Canonical **full-general UV** renormalisable one-loop RGEs from `1e7` to `1e5` GeV.
- **Remove fermions**, project to the scalar-only EFT1. Canonical **full-general EFT1** renormalisable RGEs from `1e5` to `1e3` GeV.
- Evaluate existing general Matchete **conditional one-scale** C5 Yukawa-channel kernel at the scalar threshold. B(-1)/C(-1) necessarily use the pre-existing **different conditional mixed-scalar** evaluator and its additional fixed assumptions.
- Run SM + Weinberg coefficient from `1e3` to `1e2` GeV and compute three neutrino masses and five oscillation observables, including chi2.

The **fully sequential matching** chain (all fermion EFT Wilson coefficients, their running, then all scalar matching into SMEFT) is **not implemented/certified by these existing conditional kernels**. The physical Matchete-RGBeta coupling dictionary also remains uncertified. Do not interpret cross-adapter coefficient equality as established. Block unsupported or missing models; do not replace their C5 with an approximate restricted result.

## Numerical comparisons

For each `s` and target model `i`, evaluate low-scale matrix `C5(s,i)`. When the reference `C5(s,s)` is available, calculate

`k(s,i) = Tr(C5(s,s)^dag C5(s,i)) / Tr(C5(s,s)^dag C5(s,s))`,

`epsilon(s,i) = ||C5(s,i) - k C5(s,s)||_F / ||C5(s,i)||_F`.

Also do direct **within-sweep pair comparisons** reproducing the Section 5.1 organisation: A(-2)/D(-2), B(-3)/C(-3), and B(-1)/C(-1). The `six_pair_comparison.csv` produces 2x3 rows, reporting the best complex scale and residual for each pair under both source Yukawas.

For each target `i`, compare **between** D(-2) versus A(-4) source Yukawas at the same target representation in `same_model_texture_comparison.csv`. This quantifies *benchmark dependence*, not representation dependence.

Original-source fit results are re-evaluated from the saved Yukawa matrices. If the new source-model C5 does not reproduce the saved source C5 to relative difference `<1e-5`, the source row is blocked rather than silently used as a reference.

## Installation and run (PowerShell)

From the local checkout root:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "$env:USERPROFILE\Downloads\T3_general_Dm2_Am4_dual_benchmark.zip" `
  -DestinationPath . -Force

python -X utf8 -m scripts.RunT3GeneralDualBenchmarks --self-test

python -X utf8 -m scripts.RunT3GeneralDualBenchmarks `
  --run-id thesis51_Dm2_Am4_v1 --preflight

python -u -X utf8 -m scripts.RunT3GeneralDualBenchmarks `
  --run-id thesis51_Dm2_Am4_v1

# Resume after input fixes, preserving successful model results:
python -u -X utf8 -m scripts.RunT3GeneralDualBenchmarks `
  --run-id thesis51_Dm2_Am4_v1 --resume
```

A preflight does not integrate the RGEs; it checks required saved general inputs. Missing inputs are recorded as `Blocked`. D(0) may remain blocked due to an unsupported negative self-bilinear scalar rotation; this is a physics matching limitation, not a poorly optimised Yukawa texture.

To inspect the paired patterns:

```powershell
$folder = '.\output\studies\general_scalar_extended\shared_general_dual\thesis51_Dm2_Am4_v1'
Import-Csv "$folder\six_pair_comparison.csv" | Format-Table -AutoSize
Import-Csv "$folder\summary.csv" | Format-Table source_model,model,status,chi2,epsilon_to_source -AutoSize
Import-Csv "$folder\same_model_texture_comparison.csv" | Format-Table -AutoSize
```

## Output destinations

Installed source file:
- `scripts/RunT3GeneralDualBenchmarks.py` (**full source file**, does not overwrite existing production scripts).

Generated numerical files (new run only):
- `output/studies/general_scalar_extended/shared_general_dual/thesis51_Dm2_Am4_v1/settings.json`
- `.../summary.csv` and `summary.json`
- `.../six_pair_comparison.csv` (original Section 5.1 three pairs, two benchmarks)
- `.../same_model_texture_comparison.csv`
- `.../source_independent_fit_preview.json` (pre-existing fitted source comparison, not new sweep)
- `.../D_m2/<model_tag>/{status,result}.json` and `.../A_m4/<model_tag>/{status,result}.json`
- `.../fig_dual_fixed_yukawa_chi2.png/.pdf`
- `.../fig_dual_C5_epsilon.png/.pdf`
- `.../fig_dual_theta23.png/.pdf`
- `.../fig_six_representation_pairs.png/.pdf`

The zip also contains **precomputed, source-backed preview** files under `Reports/output/thesis_general_5_1/dual_benchmark_source_preview/`. These use only the saved *independently fitted* D(-2) and A(-4) C5 matrices and **must not** be presented as the 32 new shared-input evaluations.

### Interpretation

Small `epsilon_to_source` **across both sweeps**: evidence for approximately representation-proportional C5 in these two benchmark textures, *conditional on matching assumptions*. Large discrepancies in pairwise epsilons **between sweeps**: evidence of texture-dependent interference with additional allowed Yukawa channels, or differing RGE evolution. Compare the additional-interaction inventory before attributing a difference purely to SU(2) group factors.

No independent fit is performed on target models; a high chi2 after Yukawa transfer is not evidence the target model cannot fit oscillation data with its own fitted Yukawas.
