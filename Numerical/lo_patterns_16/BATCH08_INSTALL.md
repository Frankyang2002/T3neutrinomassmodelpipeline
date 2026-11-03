# Chat A Batch 8 — Installation and local validation

All file paths are **relative to** `C:\thesistemp\T3neutrinomassmodelpipeline`. Install the ZIP at the repository root. Batch 8 adds only the following five files:

- `Numerical/lo_patterns_16/splitting_convergence.py`
- `Numerical/lo_patterns_16/BATCH08_THEORY.md`
- `Numerical/lo_patterns_16/BATCH08_INSTALL.md`
- `scripts/VerifyLOScalarSplitting16.py`
- `tests/test_lo_splitting_convergence_16.py`
 Existing Batch 1–7 code, frozen benchmark inputs and output results are not modified.

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = Join-Path $HOME 'Downloads\Chat_A_Batch08_Splitting_Convergence.zip'
Expand-Archive -LiteralPath $zip -DestinationPath $repo
Set-Location $repo

.\.venv\Scripts\python.exe -m pytest -q tests\test_lo_splitting_convergence_16.py
.\.venv\Scripts\python.exe scripts\VerifyLOScalarSplitting16.py
```

Optional check together with previous tests:

```powershell
.\.venv\Scripts\python.exe -m pytest -q tests\test_lo_kernel_factorization_16.py tests\test_lo_factorization_theorem_16.py tests\test_lo_splitting_convergence_16.py
```

New output (never overwrite):

`output/general_observable/lo_pattern_study16/analysis_batch08_01/`

- `BATCH08_SPLITTING_CONVERGENCE.md` — standalone analytical report
- `batch08_convergence.json` — all 16 × 6 model/t points and source manifest hashes
- `batch08_scaling_points.csv` — full tangent-remainder and transverse-response traces
- `batch08_model_summary.csv` — measured local power laws, endpoint nonproportionality
- `batch08_derivative_validation.csv` — saved endpoints, derivative check, finite-mixing vs leading derivative
- `SOURCE_MANIFEST.json`, `HANDOFF.json` — source-status and hash chain
- If matplotlib is installed: `batch08_quadratic_convergence.png`, `batch08_transverse_response.png`

The runner refuses a nonempty output directory. To repeat a run, supply a new directory **inside** the Chat A analysis root using `--output-dir`, for example:

```powershell
.\.venv\Scripts\python.exe scripts\VerifyLOScalarSplitting16.py --output-dir output/general_observable/lo_pattern_study16/analysis_batch08_02
```

**Validation distinction:** unit tests and a synthetic all-16 end-to-end fixture can be executed without your local frozen sources; the source-specific numerical results are not known until you run the command in your checkout. No 4-spinor physical normalization certification is asserted.
