Study B: heavy fermion mass scan

New file: Numerical/orchestration/RunHeavyFermionScaleComparison.py
Requires the previous local Study A files:
  Numerical/orchestration/RunCommonScalarComparison.py
  Numerical/orchestration/RunCommonScalarRGControls.py
and the local fitted configurations, RGE exports and matching artifacts.

PowerShell from repository root:
  python -m Numerical.orchestration.RunHeavyFermionScaleComparison
  python -m Numerical.orchestration.RunHeavyFermionScaleComparison --run --limit 2
  python -m Numerical.orchestration.RunHeavyFermionScaleComparison --run
  python -m Numerical.orchestration.RunHeavyFermionScaleComparison --run --force

Default: MF and fermion threshold multiplied together by 0.5, 1, 2;
reference Yukawa targets 0.01, 1; all eligible models.
Alternative: --mass-factors 0.5 1 1.5 --reference-yukawas 0.01

Output: Reports/output/full/heavy_fermion_scale_comparison/smallY_smallL/

Frozen control holds ALL couplings at UV values, retains one-loop matched C5.
Mass-scaling changes loop functions directly; the running-minus-frozen
comparison isolates the extra RGE effect at each mass choice. This is not
full fixed-order threshold matching. Cases violating mu_UV>mu_F>mu_S>mu_low
are marked Failed rather than silently altered.
