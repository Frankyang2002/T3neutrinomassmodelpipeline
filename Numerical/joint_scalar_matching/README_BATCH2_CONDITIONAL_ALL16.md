# Batch 2 — final all-16 *conditional* T3 neutrino-mass benchmark

**EXPLORATORY / CONDITIONAL — NOT A COMPLETE PHYSICAL ONE-LOOP NEUTRINO-MASS PREDICTION**

## Scope and matching sequence

This release **does not repeat any matching**. It reads the *existing* successful
Batch 1 corrected v6 joint-scalar output, independently validates its published
complex matrices and numerical diagnostics, and consolidates those into a
portable summary.

```text
UV --[ALL BSM heavy fermions at muF = 100000 GeV]-->
EFT1 (both BSM scalars active; all SM fermions remain active) --
[BOTH BSM scalar multiplets together at muS = 1000 GeV]--> SMEFT
```

All 16 original models are reported in source order, including ordinary,
B/C mixed, and A/D/E split-real scalar cases. No scalar is integrated out
separately. No new couplings, scalar masses, CG factors, or spinor phases are
invented.

## Install the full replacement/new files

Repository: `C:\thesistemp\T3neutrinomassmodelpipeline`

ZIP entries (all project-relative):

- `scripts/ConsolidateT3JointScalarBatch2.py` — new standalone consolidator.
- `tests/test_t3_joint_scalar_batch2_consolidation.py` — four focused synthetic contract tests.
- `Numerical/joint_scalar_matching/README_BATCH2_CONDITIONAL_ALL16.md` — this documentation.

Extract the ZIP from Downloads into the repository root, then run:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_Batch2_Conditional_All16.zip" `
  -DestinationPath . -Force

python -m unittest discover -s tests -p 'test_t3_joint_scalar_batch2_consolidation.py' -v

python -m scripts.ConsolidateT3JointScalarBatch2 `
  --input-dir 'output\general_observable\lo_joint_scalar_matching\exploratory_all16_v6_native_batch1_local' `
  --output-dir 'output\general_observable\lo_joint_scalar_matching\conditional_all16_batch2_v1'
```

The command requires existing Batch 1 `summary.json`, every valid model's
`C5_per_unknown_NGamma_exploratory/<slug>.json`, and the matching local
`numeric_backgrounds/<slug>.json` and `.npz` files. It does not trust the
absolute Windows paths recorded inside source JSON from other machines.
It reconstructs paths from the supplied `--input-dir` and model names.

**Outputs (exactly two new data files):**

```text
output/general_observable/lo_joint_scalar_matching/conditional_all16_batch2_v1/
  summary.json
  all16_conditional_mass_summary.csv
```

The runner refuses to overwrite any non-empty output directory. If you
rerun, choose `conditional_all16_batch2_v2` or another fresh directory.
The original v6 and corrected Batch 1 directories are **read-only inputs**.

## Fields in the consolidated data

The `summary.json` retains every model's status and, when available:

- Complex symmetric `3 x 3` `conditional_C5_per_N_GeV_inv` (`real`, `imag` matrices).
- Complex symmetric `3 x 3` `conditional_mnu_per_N_eV` (`real`, `imag` matrices).
- Three nonnegative descending SVD/Takagi singular values in
  `mass_singular_values_eV_per_abs_N`.
- The explicit `mS1Sq`, `mS2Sq` inputs at 1 TeV, and special scalar
  quadratic `mu2` when present, as supplied by the v6 source.
- Scalar eigenstate masses at the common 1 TeV matching threshold.
- Source-result/background hashes, original running snapshot hash, quartic
  seed hash and fermion-boundary hash, plus the candidate quartic alias map.
- Weinberg projection residual, matching-scale cancellation residual,
  SU(2) source parity and finite-value/matrix-symmetry checks.
- A by-category comparison for ordinary, mixed, self-split S1, self-split S2.

The **single CSV** includes all 16 statuses, the singular values, and
complex matrices as JSON arrays inside CSV cells. It can be reopened with
Python's standard `csv` module; `json.loads(cell)` decodes a matrix cell.

The mass convention is

\[
 m_\nu = \frac{v^2}{2} C_5,\quad v=246.22\,{\rm GeV},
 \quad 1\,{\rm GeV} = 10^9\,{\rm eV}.
\]

Everything is **per unknown signed** `N_GammaCC` for complex matrices,
and per its **absolute** magnitude for nonnegative singular values.
These numbers are *not* physical neutrino mass predictions.

## Validation and limitations

The validator checks (using the original v6 arrays, not a new loop engine):

1. All-16 input manifest, flags and model statuses.
2. Input background `.json` and `.npz` hashes against source result hashes.
3. `C5_total = C5_leading_log + C5_joint_scalar_bubble` at the correct
   ~1e-15 GeV^-1 scale (avoids NumPy's inappropriate default absolute tolerance).
4. Complex **symmetric**, not Hermitian, C5 and finite entries.
5. Original full SU(2) source parity and strict Ward residual bound `1e-8`.
6. Matching scale cancellation residual `<= 2e-9`.
7. Recomputed matrix norms and SVD/Takagi magnitudes versus published v6 metrics.
8. All scalar masses positive and original scale conventions unchanged.

Any invalid or missing model input is recorded as an **explicit failure**.
There is no zero substitution, invented term, or relaxed Ward projection.
The run returns exit code 2 for a partial result, while still writing both
summaries and the exact specific blockers.

Important physical limitations preserved: unknown absolute Matchete
GammaCC/Weyl sign and normalization; provisional Matchete-to-RGBeta quartic
aliases; illustrative running scalar mass/parameter trajectories; exploratory
muF/muS special-field identification; missing certified one-loop hard heavy-
fermion matching and EFT1 soft/counterterm subtraction; incomplete EFT1
Wilson mixing; no EWSB mass resummation. Numerical Ward residuals are
**not** an estimate of physical theoretical uncertainty.

## Validation evidence

The original corrected v6 GitHub summary was inspected and reports 16/16
conditional models and 0 Ward errors. The actual user's Windows Batch 1 run
independently reports 58 passed tests and 16/16 conditional results.

The four focused tests shipped in this ZIP passed in the preparation
environment with synthetic fixtures. **The Batch 2 consolidation has not
been run against the user's private local Windows directory here.** The
PowerShell command above performs that validation and writes the actual local
release.

Batch 3, if requested next, compares `A_-4` and `B_-1` against the existing
`Numerical/direct_one_loop/core.py` without using the full UV result as an
EFT counterterm or expecting equality of unmatched perturbative content.
