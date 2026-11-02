# Workstream 3 — final unit-reference sensitivity and numerical closeout

This is the **last numerical Workstream 3 batch** for the current exploratory
16-model benchmarks. It evaluates the precomputed 24 extra HHSS reference
responses against the **completed baseline** conditional dimension-five
matrices and writes a final closeout report. It does **not** run new EFT
thresholds, import couplings as physical, or modify prior results.

## Full files (new, install under project root)

- `Numerical/lo_certification/completeness/response_closeout.py`
- `Numerical/lo_certification/completeness/README_FINAL_SENSITIVITY_CLOSEOUT.md`
- `scripts/CloseLOPortalWorkstream.py`
- `tests/test_lo_completeness_response_closeout.py`

The previous fixed runner `scripts/RunLOPortalCompletionResponse.py` remains
unchanged. The output from its successful `portal_completion_02` run is the
input; there is no reason to rerun the expensive hard matching.

## PowerShell — local project root

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive `
  -LiteralPath "$HOME\Downloads\T3_WS3_Final_HHSS_Sensitivity_Closeout.zip" `
  -DestinationPath . -Force

python -m unittest discover -s tests -p 'test_lo_completeness_*.py' -v

# Read all saved matrices and validate the 16-model source inventory, no writes:
python -m scripts.CloseLOPortalWorkstream --check-only

# Produce ONE new immutable report; no existing source results are overwritten:
python -m scripts.CloseLOPortalWorkstream `
  --input-dir 'output\general_observable\lo_certification_completeness\portal_completion_02' `
  --output-dir 'output\general_observable\lo_certification_completeness\final_sensitivity_closeout_01'

$dir='output\general_observable\lo_certification_completeness\final_sensitivity_closeout_01'
$j=Get-Content "$dir\summary.json" -Raw | ConvertFrom-Json
$j | Select-Object status,model_count,valid_unit_reference_directions,ward_rejected_reference_directions,other_blocked_reference_directions,physical_C5_certified
Import-Csv "$dir\model_sensitivity.csv" | Format-Table -AutoSize
Import-Csv "$dir\reference_sensitivity.csv" |
  Sort-Object {[double]$_.max_relative_to_baseline_per_unit_coefficient} -Descending |
  Select-Object -First 12 model,reference_id,status,max_relative_to_baseline_per_unit_coefficient,single_reference_coefficient_for_10pct_baseline |
  Format-Table -AutoSize
Get-Content "$dir\WORKSTREAM3_CLOSEOUT.md"
```

Change `final_sensitivity_closeout_01` to a **new** run name if it exists.

## Exactly what is computed

For each extra invariant reference direction, the earlier Workstream 3 run
stored two *projected diagnostic only* flavour matrices

\[
 A_r = \frac{\partial C_5}{\partial \operatorname{Re}c_r},
 \qquad
 B_r = \frac{\partial C_5}{\partial \operatorname{Im}c_r},
\]

where the derivative uses the **declared V3 unit reference tensor** in the
current exploratory scalar/fermion convention. With hypothetical
\(|c_r|\leq\kappa\), the maximal Frobenius norm shift of this one direction is

\[
\max_{x^2+y^2\leq\kappa^2}\frac{\|x A_r+y B_r\|_F}{\|C_{5,\mathrm{baseline}}\|_F}
 = \kappa\frac{\sqrt{\lambda_\max(G_r)}}{\|C_{5,\mathrm{baseline}}\|_F},
\quad
(G_r)_{ij}=\operatorname{Re}\operatorname{Tr}(M_i^\dagger M_j),
\quad M=(A_r,B_r).
\]

This uses **both real and imaginary source responses** and accounts for their
interference, rather than adding the two matrix norms. The total listed
hypothetical upper envelope is the conservative sum of these per-reference
maxima. The default sensitivity scenarios are
\(|c_r|\leq 0.001,0.01,0.1\); change only for scenario exploration using
`--caps '0.001,0.01,0.1'`.

A coefficient threshold such as `single_reference_coefficient_for_10pct_baseline`
is the *hypothetical magnitude of the V3 reference coefficient* that could
produce a 10% Frobenius response; it is not the actual RGBeta coupling,
not an inferred bound, and not a 10% neutrino mass eigenvalue shift.

Since the relative ratios share the same unknown overall signed
`N_GammaCC` factor *within the consistent conditional basis*, that factor
formally cancels in these restricted relative comparisons. This does not
resolve the absolute normalization or physical reference coefficient phases.

## Source integrity, numerical scope, and fail-closed policy

- Requires the saved `portal_completion_02/summary.json` to report **16/16**
  baseline Ward passes, **14** added baseline directions, and **24** extra
  unit reference directions; checks model names and physical flags.
- Reads the **saved baseline 3x3 complex symmetric `dim5_sum` matrix**, not a
  hard-only matrix; checks saved norms against the stored full matrix.
- Reads **both** `Re` and `Im` 3x3 unit-reference matrices and verifies
  symmetry, finiteness, and saved Ward residuals.
- A `WARD_REJECTED_DO_NOT_USE_AS_C5` or `REFERENCE_NUMERICAL_BLOCKED`
  direction is **not** used as a valid Weinberg response and produces an
  unavailable (not zero) full-reference envelope for its model.
- No extra coefficient is assigned to the 24 named couplings. Each row retains
  the unassigned source symbols. A named registry entry is not assumed to map
  to one definite reference basis vector.
- Each model input file's SHA-256 is recorded; all 16 source model files and
  their per-model inventories must exist. The CLI reads only existing inputs
  and writes only to a **new** child of
  `output/general_observable/lo_certification_completeness/`.
- This analyzes the **existing conditional results**, not independent UV
  amplitudes and not a new full EFT1 RGE integration.

### Output files

- `summary.json`: status, all 16 model summaries, per-invariant responses,
  hypotheses, blockers and SHA-256 source provenance.
- `model_sensitivity.csv`: baseline norm, counts, worst individual unit
  response and conditional scenario envelopes.
- `reference_sensitivity.csv`: one row per extra reference direction,
  statuses and numerical sensitivity thresholds.
- `WORKSTREAM3_CLOSEOUT.md`: self-contained source-scoped scientific report.

**Physics remains conditional**: no complete physical one-loop `C5` or
neutrino-mass certification. Remaining independent blockers include signed
Matchete `GammaCC` → Weyl normalization, named RGBeta/Matchete coupling maps
and coefficient evolution, scalar basis phases in special models, other EFT1
operator/mixing terms, finite one-loop hard matching and counterterms, and
four previously failing independent legacy scalar projections.

The **EFT ordering is unchanged:** UV → integrate out **all heavy fermions
at the same \(\mu_F\)** → EFT1 with both scalar multiplets and SM fermions
→ integrate out **both scalars at the same \(\mu_S\)** → SMEFT.

### Verification

The package was tested against 73 local Workstream 3 unit/integration tests,
including a synthetic all-16 saved-output fixture and deliberately rejected
unit responses. **This does not substitute for running against the actual
Windows `portal_completion_02` data.**
