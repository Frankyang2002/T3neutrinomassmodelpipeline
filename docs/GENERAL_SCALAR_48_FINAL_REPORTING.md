# Frozen 48-study final benchmark reports (reporting only)

## Scientific scope and source

Source: `output/studies/general_scalar_extended/bulk_48/bulk48_v1/campaign_frozen.csv`

SHA-256 (CRLF bytes): `B3B9FE0EF2F3730AE49504B51CDD7AED4D4232E6812CBBFE0D16D78CD23ED3AF`

The generator validates this exact hash, all 16 model keys, all three scenario keys, the 48 model/study combinations, and the 20/6/22 counts before writing anything. **Never use an LF-reencoded copy as the frozen input.**

Output: `Reports/output/general_scalar_48_final/`.

The campaign is *conditional general-quartic*, **not** a complete general-scalar EFT or fully matched calculation. The registered additional renormalisable scalar quartics were included where native UV/EFT1 beta exports existed. General holomorphic quadratic extensions, full general-potential sequential Wilson matching, and physically certified Matchete/RGBeta normalisations are not included. Four extreme-hypercharge assignments rely on legacy C5; B(-1) and C(-1) accepted cases are imported conditional fits, not new bulk refits. The two running-UV B(-1)/C(-1) acceptances have altered topology-quartic boundary conditions. E-class values retain the half-normalisation encoded in the CSV.

## Included files

- `scripts/BuildGeneralScalar48FinalReports.py`: new full standalone benchmark reporting script. No modifications to existing physics scripts.
- `scripts/BuildBCSixSavedRunning.py`: optional additional presentation-only export from the *already saved* six-fit B(-1)/C(-1) Batch18 validation trajectory CSV. No re-integration.
- `docs/GENERAL_SCALAR_48_FINAL_REPORTING.md`: this instruction and source-accounting document.
- `Reports/output/general_scalar_48_final/index.html`: offline interactive 16x3 status matrix, case filters, exact stored observables, and endpoint scatter.
- `Reports/output/general_scalar_48_final/models/*.html`: 16 model-specific reports, each linking the existing UV/EFT1/EFT2 Lagrangians, RGEs, group-factor TeX and new candidate extra-quartic symbolic sources.
- `Reports/output/general_scalar_48_final/benchmark_atlas.pdf`: printable summary covering all 48 status/fit rows and five neutrino observables.
- `Reports/output/general_scalar_48_final/figures/*.pdf`: editable vector PDF status matrix, accepted chi-square chart and five-observable figure. PNG companions are included.
- `Reports/output/general_scalar_48_final/tables/benchmark_comparison.csv`: numerical results for all 48 rows.
- `Reports/output/general_scalar_48_final/tables/model_inventory.csv`: representation dimensions, alpha and registered extra-quartic counts.
- `Reports/output/general_scalar_48_final/tables/saved_yukawas.csv` and `yukawa_availability.csv`: entries extracted only from saved local `fitted_config.json`/`fit.json` files in the frozen bulk tree. Missing imported B/C sources are not guessed or substituted. The optional Batch18/B(+1) collector places SHA-verified/imported matrices in `tables/imported_saved_yukawas.csv` and annotates the interactive index.
- `Reports/output/general_scalar_48_final/provenance.json`: input fingerprint, numerical/reporting scope and missing-data flags.

The packaged PDF/HTML preview was built from the exact frozen CSV snapshot and repository-verified model metadata, **without access to your Windows supplementary fit directories**. It marks saved Yukawa entries and local quartic-registry verification unavailable in the preview. **Regenerate it inside the Windows local repository** to populate available fitting matrices and local analytical links; the generator performs no numerical integrations or fitting.

## Install and build in PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_GeneralScalar48_Final_Reporting.zip" -DestinationPath . -Force
python -m pip install matplotlib reportlab
python -m scripts.BuildGeneralScalar48FinalReports --verify-only
python -m scripts.BuildGeneralScalar48FinalReports
# Optional: reuse the six existing validated B(-1)/C(-1) scalar RGE trajectory samples;
# also export source-verified imported B/C and B(+1) accepted Yukawas.
python -m scripts.BuildBCSixSavedRunning
Start-Process .\Reports\output\general_scalar_48_final\index.html
Start-Process .\Reports\output\general_scalar_48_final\benchmark_atlas.pdf
Start-Process .\Reports\output\general_scalar_48_final\bc6_saved_running.html
```

Run `--verify-only` before generation. The command fails if `campaign_frozen.csv` has been changed, even if only its newline encoding changed.

Generator read-only sources: `campaign_frozen.csv`, optional model-specific saved `fit.json` / `fitted_config.json`, and `RGE/general/generated/general_rgbeta_quartic_registry.json`. The only writable path is the dedicated new report output directory (not the frozen input). It does not invoke other modules, solvers, optimisers, Mathematica, or external APIs.

## Existing historical reports retained as historical

Analytical sources:
- `Reports/output/analytical/Lagrangian/UV.tex`, `EFT_1_after_F.tex`, `EFT_2_after_S1_S2.tex`, `C5.tex` and compiled PDFs when available.
- `Reports/output/analytical/RGE/UV.tex`, `EFT_1_after_F.tex`, `EFT.tex` and PDFs.
- `Reports/output/analytical/GroupFactors/GF_UV.tex`, `GF_EFT_1_after_F.tex`, `GF_EFT_2_after_S1_S2.tex` and PDFs.
- `Reports/output/models/T3-*/{Lagrangian,RGE,GroupFactors}/*.tex` for all 16 cases.

Historical dashboard references:
- `Reports/output/display_fermion_start/5e-9/reused/optimised/interactive_rge.html`
- `Reports/output/display_fermion_start/5e-9/reused/optimised/sensitivity_ovals.html`
- `Reports/output/display_lambda_5e9/optimised_rge.html`
- `Reports/output/thesis_assets/figures/fig08_t3b_weinberg_mass_splitting_running.pdf`

These existing files **have not been recalculated** for the general-quartic campaign. Do not label them as new scalar-sector benchmark curves.

## Missing data and what would be necessary for genuine replacements

`campaign_frozen.csv` has three scale endpoints per model and no time/scale-series arrays. The bulk fit JSONs store final fit data and solver counters but not complete per-scale UV/EFT1/final-Weinberg trajectories. **Exception:** a separate independently validated B(-1)/C(-1) six-fit Batch18 output retains `trajectory_samples.csv` for its six cases. The optional `BuildBCSixSavedRunning.py` renders those existing sampled scalar masses directly to interactive HTML and vector PDF, without rerunning physics. The existing sensitivity ovals belong to older historical sources; no corresponding 48-study parameter-variation response/covariance series is saved in the frozen campaign folder.

To make legitimate new UV/EFT1/final-Weinberg **running curves**, one must separately export the numerical trajectories from the already accepted, frozen fitted configurations using the *implemented* candidate general-quartic RGE flow. This is a deterministic forward evaluation, not another fit, but it has **not** been executed by this report generator. Genuine benchmark sensitivity ovals additionally require saved parameter variations; the central points alone are insufficient.

Do not call the legacy `Numerical.plotting.BuildDisplayResults` on `output/full` expecting general-quartic trajectories: that script reads older `output/full/comparison/*/running_diagnostics.json` and will rebuild historical plots.
