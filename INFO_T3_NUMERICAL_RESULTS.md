# T3 fixed-comparison numerical results extraction

## Verified run status (2026-10-09)

User's numerical-only resume log: `63 Success`, `1 AlreadyComplete`, `0 Failed` across four comparison scenarios (16 representations/hypercharge choices each). This is completion of numerical evaluation, **not** agreement with neutrino data. Matching/RGE artifacts were reused; Matchete was not rerun during resume.

## Extractor

`Numerical/diagnostics/ExtractFullComparisonResults.py` reads existing `output/full/comparison/<scenario>/<model_key>/<class_key>/data/running_diagnostics.json` and writes:

- `Reports/output/full/t3_full_comparison_observables.csv`
- `Reports/output/full/t3_full_comparison_observables.json`

One record per successful diagnostic, with masses (eV), two mass-squared splittings plus delta m32² (eV²), three sin² mixing angles, ordering, Takagi residual, and model/scenario identifiers. It **does not** run the symbolic/numerical pipeline, modify input outputs, or claim experimental compatibility.

Angles are reconstructed from `low_energy.pmns_abs`, using the implementation convention in `Numerical/diagnostics/RunningDiagnostics.py`:

- `sin² θ13 = |U_e3|²`
- `sin² θ12 = |U_e2|²/(1-|U_e3|²)`
- `sin² θ23 = |U_μ3|²/(1-|U_e3|²)`

Masses and splittings are taken directly from `low_energy`, not re-derived. Cases missing `running_diagnostics.json` are absent from the output; malformed/failed diagnostics are recorded under JSON `errors`. Exit code 0 only when exactly 64 valid cases exist, 16 in each expected scenario, with no malformed diagnostics.

## Commands (PowerShell, repository root)

```powershell
python -m pytest -q tests/test_extract_full_comparison_results.py
python -m Numerical.diagnostics.ExtractFullComparisonResults
Import-Csv Reports\output\full\t3_full_comparison_observables.csv | Select-Object -First 5 | Format-Table
```

## Caution

The uploaded single-case diagnostic used in field verification has large neutrino masses (`~10^6 eV`) in a large-coupling benchmark. Success indicates integration and file generation, not realistic neutrino masses. Do not interpret the 64/64 run success as a fit.
