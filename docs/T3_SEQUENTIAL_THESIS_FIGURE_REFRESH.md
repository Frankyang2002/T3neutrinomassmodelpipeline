# Sequential EFT thesis figure and Appendix B refresh

## Purpose and fixed scope

The thesis uses **D(-2)** and **A(-4)** independently optimised ordinary Yukawa matrices `(Y1,Y2)` as two *shared-reference benchmarks* across all 16 models (Section 5.1 and Appendix B). This work reruns those references with the current **fermion-first** sequential soft-sector EFT calculation:

`UV -> integrate out ALL heavy fermions at muF -> EFT1(S1,S2) -> integrate out ALL scalars jointly at muS -> SM+Weinberg at muLow`.

Only the **two donor Yukawa pairs** are reoptimised. Other scalar, mass, SM and representation-dependent inputs are preserved from each saved general-scalar fitted configuration, and are not reoptimised. Model-specific `y1Conj` / `y2Conj` are set to zero **at the UV boundary** by default, matching the thesis's D(-2) Appendix B additional-Yukawa convention; `--extra-yukawas saved` is a separate, explicitly nonuniform scenario.

All actual evaluations use the native dynamic `C_LLSS` and fixed-order `LLSS -> C5` leading logarithm plus the **joint** finite scalar bubble. The old one-scale `C5` is not substituted. The hard one-loop contribution is not included; spinor normalisation and mixed-scalar basis identification remain uncertified. Every result is a **conditional** benchmark, not certified full one-loop matching.

## Installation

Extract the ZIP into the root `C:\thesistemp\T3neutrinomassmodelpipeline` without flattening the folders. It contains:

- `scripts/RunT3SequentialThesisRefresh.py`: refit exactly D(-2) and A(-4), then transfer those fixed matrices and recompute three running-stage contrasts in 16 models.
- `scripts/PlotT3SequentialThesisRefresh.py`: regenerate PDF+PNG heatmaps, appendix CSV/LaTeX tables, and proportionality/fit diagnostic plots **from completed local output only**.
- `docs/T3_SEQUENTIAL_THESIS_FIGURE_REFRESH.md`: these instructions.

Required existing files/scripts: `scripts/RunT3All16DynamicSequential.py`, `scripts/RunT3All16DynamicSequential16.py`, `scripts/CompareT3SequentialRunningFrozen16.py`, and the same native Matchete, RGBeta and saved-fit inputs used for the successful `batch16_v1` campaign.

## PowerShell commands

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath ".\T3_Sequential_Thesis_Figure_Refresh.zip" -DestinationPath . -Force

# Input/source check (no numerical solve, no writes):
python -X utf8 -m scripts.RunT3SequentialThesisRefresh --phase fit --run-id thesis_seq_v1 --check-only

# 1. Fit D(-2) and A(-4) using the new sequential C5, not the old one-scale C5:
python -X utf8 -m scripts.RunT3SequentialThesisRefresh --phase fit --run-id thesis_seq_v1 --max-nfev 60

# Review both chi2 values in the generated fits before proceeding:
Get-ChildItem .\output\general_sequential_pilot\thesis_sequential_refresh\thesis_seq_v1\fits\*.json |
    ForEach-Object { $j = Get-Content $_.FullName -Raw | ConvertFrom-Json;
      [pscustomobject]@{ Model=$j.model; Status=$j.status; Chi2=$j.final_chi2 } } |
    Format-Table -AutoSize

# 2. Transfer each fitted Y1,Y2 unchanged to all 16 models, and do UV/EFT1/SMEFT frozen cases:
python -X utf8 -m scripts.RunT3SequentialThesisRefresh --phase transfer --run-id thesis_seq_v1

# 3. Generate figures and tables; intentionally NOT sensitivity Figure 5.4:
python -X utf8 -m scripts.PlotT3SequentialThesisRefresh --run-id thesis_seq_v1
```

If a donor does **not** satisfy the project's approximate-screening `chi2 <= 5.86`, the transfer stage stops rather than presenting it as an accepted result. Adjust `--max-nfev` and use a **new** run ID to attempt a different fit (do not overwrite a completed fit). To explore a poor fit anyway, use `--include-poor` with both the transfer and plot commands; this produces provisional, explicitly labelled outputs.

The two fit steps may involve many individual RGE and matching evaluations. They run locally and do not require Matchete or GitHub writes.

## Generated inputs and figures

Physics outputs (new folder, preserves original batches):

```
output/general_sequential_pilot/thesis_sequential_refresh/thesis_seq_v1/
  manifest.json
  fits/D_m2.json
  fits/A_m4.json
  transfers/D_m2/<16 model slugs>/transfer.json
  transfers/A_m4/<16 model slugs>/transfer.json
  transfer_observables.csv
  running_effects.csv
  transfer_phase_summary.json
```

Thesis graphics and tables (new location, **never** overwrites historical PDFs):

```
Reports/output/thesis_sequential_refresh/thesis_seq_v1/
  figures/fig5_2_sequential_D_m2.pdf + .png
  figures/fig5_2_sequential_A_m4.pdf + .png
  figures/fig5_3_sequential_D_m2.pdf + .png
  figures/fig5_3_sequential_A_m4.pdf + .png
  figures/fig_sequential_c5_proportionality.pdf + .png
  figures/fig_sequential_two_donor_fits.pdf + .png
  tables/table_B1_Dm2_observables.csv + .tex
  tables/table_B2_Am4_observables.csv + .tex
  tables/table_B3_pmns_shared.csv + .tex
  tables/table_B4_c5_proportionality.csv + .tex
  captions_and_scope.md
  figure_manifest.json
```

**Figure 5.1** is the analytical B/C additional one-loop diagram: unchanged, no new Yukawa-dependent data needed. **Figure 5.4** is the coupling sensitivity graph: explicitly excluded and untouched. Group-theory Tables 5.1 and 5.2 are also analytical and unchanged.

## One-stage-off definition

Every case matches at `muF = 10^5 GeV` and `muS = 10^3 GeV`. The `C_LLSS -> C5` fixed-order logarithm remains ON in all cases. The main full-running evaluation is `U1_E1_S1`.

- UV effect: `100*(X(U1_E1_S1) - X(U0_E1_S1))/X(U0_E1_S1)`.
- EFT1 effect: `100*(X(U1_E1_S1) - X(U1_E0_S1))/X(U1_E0_S1)`.
- SMEFT effect: `100*(X(U1_E1_S1) - X(U1_E1_S0))/X(U1_E1_S0)`.

`U0_E0_S1` is also stored to support the four-case comparison, but is not one of the three heatmap contrasts.

The historical thesis Section 5.2 heatmaps used 16 separately optimised Yukawa points. These new heatmaps deliberately use the **same D(-2) or A(-4) Yukawa pair in every model**; they should be described as *shared-reference running heatmaps* in Chapter 5, not a numerically exact reproduction of the old individual-fit figures.

In the current mixed B(-1)/C(-1) implementation, the `U(muS)` scalar rotation is used to identify the F-boundary basis. This remains exploratory even if numerical scale-cancellation residuals are small.

## No uncomputed or invented figures

The packaged Python files generate the figures **only after successful local fits and 32 transfer calculations**. The package does not contain fabricated or refitted benchmark results. It was syntax-checked and exercised with synthetic input data to verify the plotting and output formats; those test figures are *not* scientific outputs and are **not included** in the ZIP.
