# Batch 2: read-only 16-model numerical interpretation

This file accompanies `Numerical/lo_patterns_16/result_analysis.py`. It consumes
only the **existing** `provisional_01/` output of `RunLOPatternStudy16.py` and
writes an entirely new `analysis_batch02_01/` directory. It does not recompute
UV masses or alter any of the production pipeline, benchmark inputs, or prior results.

## What is computed

- A table of all **three** sorted Takagi singular values per model, plus the
  sum and Frobenius norm of the full complex symmetric neutrino mass matrix.
- Actual frozen-benchmark family means, ranges, and hypercharge-dependence.
- Exact polynomial checks for scans of one Yukawa at a time. Write
  `M11 = mnu[y1_HC x y1_HC]`, `M12 = mnu[y1_HC x y2_HC]`, and
  `M22 = mnu[y2_HC x y2_HC]`. For fixed masses and scalar mixing,
  
  \[
  m_\nu(ty_1,y_2)=t^2 M_{11}+t M_{12}+M_{22},\qquad
  m_\nu(y_1,ty_2)=M_{11}+t M_{12}+t^2 M_{22}.
  \]

  The script compares the singular values predicted by these exact identities
  against the **persisted** scan rows. Simple linear Yukawa scaling is additionally
  reported, but only holds when the same-vertex terms vanish.
- Scalar-mixing sign-flip, zero limit and `t=2` linear-approximation residuals.
  For the restricted `B=0` two-scalar-block topology the complex mass is odd in
  the off-diagonal mixing and the Takagi masses are even under its sign flip.
  Finite mixing need not scale exactly linearly.
- Heavy-mass scan: empirical log slope between factors 0.5 and 2; it is **not**
  forced to equal `-1` because the full loop functions depend on squared-mass ratios.
- Per-charge matrix norms and, for models with multiple channels, real and
  imaginary Frobenius overlaps. The coherence ratio
  \[\mathcal C=\frac{\|\sum_Q M_Q\|_F}{\sum_Q\|M_Q\|_F}\]
  is at most 1 (up to roundoff), and diagnoses constructive/destructive matrix
  interference. The report includes a remainder term if the source-CG topology
  charges do not exhaust the complete component charge decomposition.

**Do not add Takagi eigenmasses channel by channel.** The full mass matrix is
additive, its singular values are not. The group-only weights from
`common_kernel_group_only.json` remain a separate artificial common-kernel
comparison, not a model-independent physical mass prediction.

## PowerShell commands

From `C:\thesistemp\T3neutrinomassmodelpipeline` after extracting the Batch-2 ZIP:

```powershell
.\.venv\Scripts\python.exe -m pytest -q tests\test_lo_pattern_analysis_16.py

.\.venv\Scripts\python.exe scripts\AnalyseLOPatternStudy16.py
```

Inputs (unchanged):
`output\general_observable\lo_pattern_study16\provisional_01\`.

Outputs (must be new/empty):
`output\general_observable\lo_pattern_study16\analysis_batch02_01\`.

For a repeat analysis use a new directory, for example:

```powershell
.\.venv\Scripts\python.exe scripts\AnalyseLOPatternStudy16.py `
  --output-dir 'output\general_observable\lo_pattern_study16\analysis_batch02_02'
```

`SOURCE_MANIFEST.json` captures SHA256 of the prior source results. All frozen
benchmark, provisional component-input, source-audit and executable guard
SHA256s are checked against the prior run before any reports are written.
A mismatch, blocked/missing scan, incomplete model inventory, false physical
certificate or inconsistent matrix/scan identity fails closed.

## Expected output artifacts

- `BATCH02_REPORT.md`: complete 16-row comparison and conditional interpretation.
- `batch02_analysis.json`: machine-readable metrics and within-family patterns.
- `batch02_model_summary.csv`: all 16 rows with masses and scaling diagnostics.
- `batch02_charge_matrix_norms.csv`: charge-pair partition contributions.
- `batch02_numerical_checks.json`: exact-polynomial checks.
- `SOURCE_MANIFEST.json`: prior result SHA256s.
- `HANDOFF.md`, `HANDOFF.json`: batch provenance and pending B1/B2 evidence.
- Up to four PNG plots if matplotlib is installed.

**Limits:** This batch cannot upgrade 0/16 physical normalization certification.
The numerical tests use a fully synthetic report schema. They do not replace
16 actual frozen-benchmark evaluations, which the previous run already completed.
No EFT work, running or threshold matching is performed.
