# Fermion-start reports via established plot builders

This study does NOT replace `InteractiveModelComparison.py` or `BuildFittedModelEnvelopes.py`.
The adapter `Numerical/orchestration/BuildFermionStartExistingReports.py` prepares saved data in their input structures and then invokes these exact pre-existing HTML renderers.

## Numerical convention

`011`: input T3 UV-state couplings are specified at the fermion threshold mu_F;
UV evolution is frozen, fermion matching enabled, intermediate scalar EFT RGEs enabled,
scalar matching enabled, and final SM+Weinberg RGEs enabled.
**The 'UV' interval is not displayed as a physical running curve.**

Cases `5e-9`, `5e-10`; class E has lambdaT3 halved in the fitted configuration.
Modes `optimised`, `shared_a_m2`. Neither fitting nor symbolic RGEs are rerun.

## Run

```powershell
# Verify a single fitted model and established HTML builders first
python -m Numerical.orchestration.BuildFermionStartExistingReports `
  --case 5e-9 --mode optimised --model T3-A-m2 --kind all

# Generate all cases, models, and modes (about 2*2*16*(3*11+1) numerical evaluations)
python -m Numerical.orchestration.BuildFermionStartExistingReports `
  --case all --mode all --kind all
```

Files are saved under:

- `output/studies/fermion_start_lambda_comparison/<case>/reports_existing/<mode>/models/<model>/running_diagnostics.json`
- `output/studies/fermion_start_lambda_comparison/<case>/reports_existing/<mode>/sensitivity/<model>/sensitivity_diagnostics.json`
- `Reports/output/display_fermion_start/<case>/reused/<mode>/interactive_rge.html`
- `Reports/output/display_fermion_start/<case>/reused/<mode>/sensitivity_ovals.html`

The existing interactive RGE renderer consumes saved final-Weinberg curves. Its UV quantities are deliberately unavailable and will not plot. Scalar-only EFT1 state-vector samples and their exact block layout are saved for later plotting, without manufacturing an effective Weinberg coefficient above mu_S.

The existing fitted envelope renderer uses saved gamma directions T1/T2/T3. Their definitions are preserved from `Numerical/fitting/BenchmarkSensitivity.py`, but each evaluation uses the fermion-start 011 physics prescription instead of the previous full-UV evaluator. Ovals are empirical, descriptive envelopes, **not confidence contours**. Its NuFIT/JUNO overlays retain the renderer's existing definitions.

Default sensitivity gamma range is 0.1 to 10 with 11 log points; change to e.g. `--gamma-min .001 --gamma-max 1000 --points 25 --force` if the wider original scan is desired. The default is a practical pilot, **not identical** to the former fitted-sensitivity gamma range. `--force` recomputes all selected saved data.

### Checks / limits

- Each RGE curve's low-energy endpoint is checked against the corresponding `comparison` evaluation.
- A gamma=1 centre is checked against the saved benchmark for every sensitivity direction.
- Failed individual gamma points are retained as failed, not interpolated.
- Existing plotting templates are reused unchanged.
- The full `InteractiveModelComparison` legacy intermediate *direct C5 continuation* is omitted because the fermion-start evaluator does not save that separate matching diagnostic. Its other available plotting selectors still work, and no curve is invented.
- The script is syntax-checked, but not physics-executed against the user's local symbolic artifacts. Run the one-model command before the all-model command.
