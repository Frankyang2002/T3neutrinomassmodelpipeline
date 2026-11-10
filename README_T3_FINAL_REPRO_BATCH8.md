# T3 final all-16 numerical reproducibility package — Batch 8

<!-- T3_CURRENT_SCOPE_BEGIN: 2026-11-10 -->

**Document scope:** Batch 8 is a historical conditional report, not current full-general output.

**Current implementation status (10 November 2026):** The requested A(-4)/D(-2) **full-general, complete-one-loop, fermion-first -> joint-scalar -> SMEFT Yukawa fit is not yet an executable end-to-end numerical workflow**. The repository already implements general Matchete/RGBeta components and an actual symbolic `RunSequentialT3Matching` stage engine, but its general Wilson matching outputs are not fully connected to the two-model numerical fit. The current thesis Sobol pilot is **restricted/conditional** (analytic fermion triangle and selected scalar bubble), even if it reaches a low chi-square. See [implementation status](docs/IMPLEMENTATION_STATUS.md), [A(-4)/D(-2) feasibility](docs/T3_A4_D2_GENERAL_PIPELINE_FEASIBILITY_2026_11_10.md) and [implementation batches](docs/T3_A4_D2_FULL_GENERAL_IMPLEMENTATION_BATCHES.md). Earlier numerical results and verification claims below refer to their **dated, historical configurations** unless an individual run manifest establishes otherwise.

<!-- T3_CURRENT_SCOPE_END -->
**Scope:** One command to resume the existing all-16 Batch-7 numerical pipeline from the saved Batches 1/2/4/5, validate that its result is reproducible, and consolidate *every* numerical branch into a stable final export. **No changes** to the existing UV/EFT1/scalar/SMEFT physics. **No changes to GitHub.**

Prepared using the locally available **Batch 7 v2 resume hotfix** and the public code at GitHub commit `4fdae67ce2ee5f22d64284a2c83ebc4ffb3400f0`, inspected 10 November 2026. Source hashes are taken from the installed local files, not assumed equal to the public commit.

## Install and run (PowerShell)

From your project directory:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "$HOME\Downloads\t3_final_repro_batch8_all16.zip" `
  -DestinationPath . -Force

# Final one-command execution: tests, Batch 7 --resume, validation, export
.\scripts\RunT3FinalRepro.ps1
```

The complete Batch-7 run already exists (16/16 in your reported output). The default script uses `--resume`, so if source hashes and numerical scales match it should select **no models for recalculation**. Existing Batch-7 results are not deleted.

The PowerShell wrapper invokes `.venv\Scripts\python.exe`. Equivalent direct invocation:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3FinalRepro.py
```

To **only reassemble from existing results** without invoking the Batch-7 runner:

```powershell
.\scripts\RunT3FinalRepro.ps1 -CollectOnly
```

To **verify the completed final report** against the current saved inputs, numerical outputs and installed code, without modifying any files:

```powershell
.\scripts\RunT3FinalRepro.ps1 -VerifyOnly
```

For independent testing of the report implementation:

```powershell
.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_final_repro_batch8.py
```

## Intended project paths (all complete new files)

| Project-relative path | Purpose |
|---|---|
| `Numerical/parallel_t3/wilson/final_repro_batch8.py` | Collection, source hash validation, C5 algebra and low-energy Takagi consistency, deterministic CSV/NPZ/JSON writing and verification |
| `scripts/RunT3FinalRepro.py` | One-command entry point; tests + Batch7 resume + consolidation |
| `scripts/RunT3FinalRepro.ps1` | Windows PowerShell launcher |
| `tests/test_t3_final_repro_batch8.py` | Isolated 16-model fixture and failure-mode tests |
| `README_T3_FINAL_REPRO_BATCH8.md` | This README |

Nothing in this package replaces a Batch 1–7 physics source file.

## Output directory

`output\parallel_t3\final_repro_all16_batch8_v1\`

| Output | Contents |
|---|---|
| `all16_overview.csv` | 16 nominal analytic-triangle models with masses, mass-square differences, mixing angles (where defined), counterterm and pole diagnostics |
| `all16_nominal.csv` | 48 rows: 16 models × 3 conditional prescriptions, at mu = MF |
| `all16_scenarios.csv` | 144 rows: 16 models × 3 prescriptions × 3 matching subtraction scales |
| `all16_wilson_matrices.npz` | Complex 3×3 C5 at MS and low scale for all 144 branches, plus per-model hard/running/scalar/CT/native totals |
| `manifest.json` | SHA256 of 14 prerequisite files per model, the Batch-7 ledger and 3 result files per model, installed source-code hashes, output hashes and numerical conventions |
| `PROCESS.md` | Reproducible EFT order, assumptions and prescription labels |

All results are **read from saved numerical tensors and rechecked**. The report does *not* fabricate or rerun Matchete/RGBeta symbolic inputs. The extraction also checks that 3 heavy fermions are removed together at MF = 100000 GeV, EFT1 contains both scalar multiplets, both scalars are removed together at MS = 1000 GeV, and SMEFT is evolved to the recorded low scale.

At the nominal subtraction scale the primary branch must satisfy

\[
C_5^{\mathrm{SMEFT}}(M_S) =
 C_5^{F\text{-hard}}(M_F)
 + \Delta C_5^{\mathrm{EFT1},\,\text{strict}}
 + C_5^{\mathrm{scalar\ bubble}}(M_S),
\]

and all three nominal scenarios must reproduce their saved Batches 2/4/5 C5 arrays. The masses are independently recalculated from the stored low-scale C5 through

\[
m_\nu=-\frac{v^2}{2}C_5, \qquad
m_i=10^9\,\sigma_i(m_\nu)\ \mathrm{eV},
\]

where the singular values are sorted ascending. This does **not** assign normal/inverted physical neutrino ordering.

## Interpretation

This package completes the **numerical engineering and reproducibility** objective. It retains the implemented conditional assumptions of Batches 1–7. `native_formal` uses the formal native pole-dropped prescription; `native_as_UV_minus_soft_hypothesis` is a diagnostic alternative, *not* established full UV minus soft matching. Nonzero native-pole residuals remain in the data. They are not silently set to zero or absorbed into artificial counterterms. The 16 benchmark masses are not physical neutrino-mass fits.

This report has been validated against synthetic saved-output fixtures, not the unshared Windows Batch-7 files. Your Windows execution above confirms Batch 7 16/16; Batch 8 still requires one final local run to confirm all saved-data layout and provenance checks.

## Failure behaviour

If a needed source file has changed, a Batch-7 output is missing or its SHA256 differs from the saved ledger, a 3×3 complex C5 is non-symmetric, a nominal matching identity fails, or low-scale masses disagree with stored C5, final generation returns nonzero and does not produce a false 16/16 success. If the report is run with `--collect-only`, prior batches are **never rerun**. All final outputs are confined to the new Batch-8 folder, and paths inside any existing Batch-1/2/4/5/7 directory are rejected for output.
