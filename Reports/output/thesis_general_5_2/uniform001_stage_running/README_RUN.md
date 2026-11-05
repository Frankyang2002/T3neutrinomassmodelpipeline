# Thesis Section 5.2 — uniform-0.01 running-stage figures

## 2026-11-05 repair: line-ending portable accepted-source checks

Some saved reference source fingerprints were calculated from Windows **CRLF**
files; later checkouts may use **LF** for the exact same JSON/WL/InputForm content.
This affected the B(-1)/C(-1) reference validator after the A and B(-3)
reference checks had passed. We verified that all mismatched B(-1)/C(-1)
source-file hashes on the GitHub main branch match the saved hashes after
**only** LF-to-CRLF conversion (including UV/EFT1 beta JSON). On a user's
local checkout the script performs the same verification independently.

The stage-running driver now tries the unchanged strict campaign validator
first. **Only** when it reports `Original accepted-source fingerprints changed`,
it checks the original path inventory, exact accepted manifest fingerprint,
and each saved SHA-256 using its unchanged bytes and its LF/CRLF equivalents.
It does not normalise numerical values, whitespace, JSON formatting, or
symbolic expressions. Any substantive file-content change still causes an
error naming the affected source key. The existing new-benchmark scalar/Yukawa
checks and full-running low-energy C5/observable reproduction remain mandatory.
No source files or previously computed numerical results are rewritten.

Install replacement files from the repo root and resume the same run ID:

```powershell
Expand-Archive -LiteralPath .\t3_uniform001_running_eol_fix.zip -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests/python -p 'test_t3_uniform001_stage_running.py' -v
python -X utf8 -m scripts.RunT3Uniform001StageRunning --stage preflight --models B_-1,C_-1
python -X utf8 -m scripts.RunT3Uniform001StageRunning --stage compute --models all --run-id uniform001_running_v1 --resume
python -X utf8 -m scripts.RunT3Uniform001StageRunning --stage plot --run-id uniform001_running_v1
```

**Do not plot** until the compute stage reports 16/16 completed and passes
all reference baselines. This remains conditional one-scale matching, not a
certified all-fermions -> all-scalars -> SMEFT calculation.


This is a **new, independent study** under the already accepted 16-model uniform-0.01 benchmark. It **does not reuse** the old `Reports/output/thesis_qualitative/representation_running/csv/frozen_stage_observable_changes.csv` rows or the old thesis heatmaps.

- The benchmark is all allowed **BSM quartics other than** `lambdaT3` equal to **+0.01 at 10^7 GeV**, allowed conjugate Yukawa matrices **zero at the UV boundary**, and 18 ordinary Yukawa entries taken from each model's **new accepted fit**, without refitting.
- Baseline `lambdaT3` is **0.1** for A–D and **0.05** for E. SM Higgs quartic, gauge/SM Yukawas, masses, scalar bilinear prescription and matching scales are unchanged.
- The UV, EFT1 and SMEFT comparisons freeze one stage at a time. The other stages use the same existing numerical solver and conditional matching as the reference fits.
- Four configurations per model: `111` (fully running); `011` (UV off); `101` (scalar-EFT running off); `110` (SM+Weinberg running off). The percentage reported is `100*(O_111-O_off)/O_off`.
- One-scale **conditional** full-general Matchete C5 is recomputed at each scalar endpoint. The calculation is **not** a completed/validated sequential EFT Wilson matching, and in particular it does **not** independently transport EFT1 LLSS Wilson coefficients. No scientific conclusion about a certified full sequential EFT should be based on these plots.
- For the 14 nonmixed models the driver uses the explicit pinned-canonical cache and baseline check from `BuildT3Uniform001FixedYukawaTables.py`. It does **not** edit live sources or claim the old cache hashes match today's Lagrangian/RGBeta sources. For B(-1)/C(-1) it uses the existing conditional mixed-scalar backend with **absolute 0.01** extra scalar quartics.
- Every computed `111` observable and low-scale C5 is compared against the saved new reference result. A mismatch blocks that model. Only 16/16 verified results can generate thesis heatmaps.

## Project-relative files in this package

- `scripts/RunT3Uniform001StageRunning.py`: complete new driver, including physics evaluation and report assembly.
- `tests/python/test_t3_uniform001_stage_running.py`: isolated validation tests.
- `Reports/output/thesis_general_5_2/uniform001_stage_running/README_RUN.md`: this document.
- `Reports/output/thesis_general_5_2/uniform001_stage_running/THESIS_FIGURES.tex`: LaTeX captions and image paths (use after plot generation).
- `Reports/output/thesis_general_5_2/uniform001_stage_running/THESIS_CONSISTENCY.md`: thesis text corrections.

## PowerShell commands (repository root)

```powershell
Expand-Archive -LiteralPath .\t3_uniform001_running_update.zip -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests/python -p "test_t3_uniform001_stage_running.py" -v

python -X utf8 -m scripts.RunT3Uniform001StageRunning --stage preflight --models pilot
python -X utf8 -m scripts.RunT3Uniform001StageRunning --stage compute --models pilot --run-id uniform001_running_v1

# Resume from the two pilot models; expensive calculations are checkpointed by model.
python -X utf8 -m scripts.RunT3Uniform001StageRunning --stage compute --models all --run-id uniform001_running_v1 --resume

# This command refuses to create either heatmap unless all 16 models are complete.
python -X utf8 -m scripts.RunT3Uniform001StageRunning --stage plot --run-id uniform001_running_v1
```

The numerical checkpoint directory will be:

`output/studies/thesis_qualitative/uniform_scalar_001/stage_running/uniform001_running_v1/`

The updated **original-style** plots will be:

- `Reports/output/thesis_general_5_2/uniform001_stage_running/uniform001_running_v1/plots/oscillation_mixing_heatmap_percent.png`
- `Reports/output/thesis_general_5_2/uniform001_stage_running/uniform001_running_v1/plots/oscillation_masssq_heatmap_percent.png`

The complete comparison CSV is `.../stage_running/uniform001_running_v1/frozen_stage_observable_changes.csv`, with 16×3×5 = **240** model-stage-observable rows.

Do not interpret `Completed` numerical evaluations as evidence of acceptable fits for the **stage-off** cases; the baseline new-reference fits are accepted, while stage-off predictions are intentionally not reoptimised.

## Important thesis consistency issue

The existing PDF's Section 3.3 presents a **sequential fermion -> scalar matching** narrative and EFT1 higher-dimensional Wilson running as if it were a completed validated numerical calculation. The current general-model scalar-0.01 results use the **conditional one-scale** coefficient. Change those implementation claims; a physically complete sequential EFT must be independently generated and validated rather than inferred from the present running-stage plots.
