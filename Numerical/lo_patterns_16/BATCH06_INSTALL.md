# Chat A Batch 6 — Install and run

## Full new project-relative files

- `Numerical/lo_patterns_16/kernel_factorization.py` — all-16 Yukawa-independent heavy coefficient contraction and comparisons
- `Numerical/lo_patterns_16/BATCH06_THEORY.md` — source-derived derivation, criteria, and limitations
- `Numerical/lo_patterns_16/BATCH06_INSTALL.md` — this file
- `scripts/ProveLOKernelFactorisation16.py` — run wrapper
- `tests/test_lo_kernel_factorization_16.py` — algebra, negative controls, synthetic 16-model test, optional local integration

**Prerequisites:** Previous Chat A Batches 1–5 installed and successfully run; unmodified candidate Weyl source inputs, frozen benchmark JSONs, and the saved `provisional_01`, `analysis_batch03_01`, `analysis_batch04_01`, and `analysis_batch05_01` outputs. Python `.venv` containing NumPy and pytest; matplotlib only for an optional plot. This package does not change GitHub, source UV modules, or earlier outputs.

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = Join-Path $HOME 'Downloads\Chat_A_Batch06_Heavy_Kernel_Factorization.zip'
Expand-Archive -LiteralPath $zip -DestinationPath $repo
Set-Location $repo
.\.venv\Scripts\python.exe -m pytest -q tests\test_lo_kernel_factorization_16.py
.\.venv\Scripts\python.exe scripts\ProveLOKernelFactorisation16.py
```

If a rerun is needed, **do not overwrite** the first output. Supply a fresh target:

```powershell
.\.venv\Scripts\python.exe scripts\ProveLOKernelFactorisation16.py `
  --output-dir output\general_observable\lo_pattern_study16\analysis_batch06_02
```

Output directory: `output/general_observable/lo_pattern_study16/analysis_batch06_01/`. Primary report: `BATCH06_ANALYTICAL_FACTORISATION.md`; numeric data: `batch06_factorization.json`; companion CSVs: `batch06_coefficient_tests.csv`, `batch06_symmetry_pairs.csv`, `batch06_flavour_matrix_tests.csv`, `batch06_linear_approximation.csv`, `batch06_E_charge_channels.csv`, and `batch06_tensor_reconstruction.csv`; provenance `SOURCE_MANIFEST.json`, handoff `HANDOFF.json`, and optional chart `batch06_factorization_stress.png`.

The normal suite includes one **skipped** opt-in test requiring the actual Windows checkout. To run the separate integration test (only if its unique target is unoccupied):

```powershell
$env:RUN_LO_KERNEL_INTEGRATION='1'
.\.venv\Scripts\python.exe -m pytest -q tests\test_lo_kernel_factorization_16.py
Remove-Item Env:\RUN_LO_KERNEL_INTEGRATION
```

Important: the coefficient-level factorization criterion is exact as an algebraic implication, but a small **numerical** residual is evidence at the explicitly sampled spectra, not a general theorem for the full UV model. Counterfactual mass changes are not validated Lagrangian couplings. Physical normalization certification remains **0/16** in this workstream.
