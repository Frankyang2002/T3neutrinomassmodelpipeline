# Thesis-to-repository artifact map

**Basis:** uploaded November 2026 109-page thesis and GitHub `main` commit `91f8d7b0f81417df751632bacb436156a04872de`. Paths below were checked for existence in that commit unless described as missing. These are locations of **possible or recorded provenance**; they are not all demonstrated to be the exact inputs of the final uploaded PDF. Do not silently replace a figure's data with a newer campaign.

| Thesis content | Keep these known code or data locations | Status and caveat |
| --- | --- | --- |
| §3 representations and Appendix C interaction inventory | `common/T3Model.py`, `Lagrangian/model/`, `Lagrangian/interactions/`, `RGE/group_factors/`, `input/lo_direct_one_loop/` | Core scientific source; no automatic removal |
| §3 UV and intermediate EFT running | `RGE/running/`, `RGE/general/`, `Numerical/core/`, `Numerical/running/`, `Numerical/parallel_t3/` | Some selected routes are conditional; preserve even if filename says `batch` |
| §3 neutrino mass and mixing | `RGE/matching/`, `physics/`, `Numerical/diagnostics/`, `Numerical/fitting/` | Scientific source and numerical conventions |
| §4 fitting and numerical benchmarks | `configs/nufit_6_1_ic24_no.json`, `Numerical/benchmarks/fitted_reference/`, `output/studies/thesis_qualitative/uniform_scalar_001/campaigns/uniform001_all16_reference_v1/reference/` | Distinct fit generations; preserve source fingerprints and Yukawa matrices |
| §5.1 and Appendix B, A(-4)/D(-2) fixed-Yukawa transfer | `scripts/BuildT3Uniform001FixedYukawaTables.py`, `Reports/output/thesis_general_5_1/uniform001_fixed_yukawa_tables/uniform001_Dm2_Am4_fixed_v1/` | Contains 32-point transfer CSV/TeX and provenance for **one** accepted benchmark generation |
| §5.2 running-stage heatmaps | `scripts/RunT3Uniform001StageRunning.py`, `Reports/output/thesis_general_5_2/uniform001_stage_running/` | On inspected main only `README_RUN.md`, `THESIS_CONSISTENCY.md`, `THESIS_FIGURES.tex` were present; regenerated `uniform001_running_v1` PNGs were **not** tracked |
| §5.2 older saved running comparisons | `Reports/output/thesis_qualitative/representation_running/`, including `csv/frozen_stage_observable_changes.csv`, `plots/oscillation_mixing_heatmap_percent.png`, `plots/oscillation_masssq_heatmap_percent.png` | Existing figure alternatives; **do not archive until final-PDF provenance confirmed** |
| §5.3 T1/T2/T3 sensitivity | `Reports/output/thesis_general_5_3/solar_atmospheric_sensitivity_v1/`, `Reports/output/thesis_general_5_3/zee_uniform001_ovals_v1/full_sensitivity/` | Both old/new campaign outputs exist; their provenance identifies distinct source data |
| Figs. 5.2–5.4 and other installed final graphics | `Result Graphs/`, `Reports/figure/`, `Numerical/plotting/`, `Numerical/analysis/BuildT3ThesisAssets.py` | Keep graphics, plot code and any referenced input manifest |
| Appendix C.1 fitted matrices | `Numerical/benchmarks/fitted_reference/`, `output/studies/thesis_qualitative/`, `configs/` | Do not substitute an unrelated campaign's Yukawas based on matching model names |

## Evidence and preservation rules

1. The user's **thesis figures, numerical tables and accepted fitted-Yukawa configurations** are primary retained artifacts. Their source JSON/CSV/WL exports, reference hash manifests, and generating Python/Wolfram scripts must remain available.
2. Preserve all scientific modules in `Lagrangian/`, `RGE/`, `Numerical/`, `physics/`, `common/`, `model/`, and their active tests while the calculation dependency closure is unresolved. A static no-import check cannot prove a standalone script was never used.
3. `output/` and `Reports/output/` are **not scratch directories**. Some values are regenerated from them; others are the only surviving provenance. Do not bulk-delete.
4. Safe initial archival candidates: unreferenced root notes (`README_BATCH*`, `WORKER_*`, etc.), historical `docs/research/` notes, and identifiable LaTeX build auxiliaries (`.aux`, `.toc`, `.out`, `.compile.log`, `.log` with a retained `.tex`). A live reference from retained tracked files blocks archival.
5. Historical studies and one-time programs can be moved to an **external** archive only after checking source references, figure provenance, CLI usage and any saved-data dependencies. A restoration manifest with hashes is required.

## Reproducibility check after any deletion

From the local checkout:

```powershell
python -X utf8 pipeline.py --help
python -X utf8 full16_symbolic.py --help
.\scripts\RunCoreTests.ps1
```

For each thesis figure/table, additionally rerun its *specific* report generator in saved-data mode if the documentation and local artifacts establish a reproducible command. These checks distinguish an intact software interface from a successful live Wolfram run or a physically complete matching computation.
