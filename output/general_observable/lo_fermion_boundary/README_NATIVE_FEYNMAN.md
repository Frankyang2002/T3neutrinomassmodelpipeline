# Native Matchete Feynman-rule evidence (A_-4)

All four files are NEW inside `output/general_observable/lo_fermion_boundary/`.
No production matching, RGE, UV Lagrangian, GitHub or certified result files are modified.

## What this actually does

- Cryptographically validate the existing A_-4 source against the unchanged
  16-model source-resolved boundary (Python preflight).
- Load the user's installed Matchete v0.6.0 via `UsingFrontEnd`.
- Rebuild the **full-general** A_-4 UV model with **F heavy and BOTH scalars light**.
- Extract only the `y1` and `y2` Yukawa terms and their conjugates, while
  validating that all general interactions were assembled.
- Apply native Matchete `FeynmanRules` separately to the UV Yukawa vertices
  and the originally matched holomorphic PL LLSS contact operator.
- Export their exact native `InputForm` when small, or a SHA and length plus
  prefix when too large; report source term and mass term separately.
- Validate the provenance, structure, and explicit non-certification of the
  JSON evidence with Python.

**Not done yet:** the two UV vertices are not contracted with a signed Dirac
propagator and both crossing/orderings. This is native vertex extraction, NOT
an independent connected four-point amplitude, and it cannot certify the
GammaCC-to-master Weyl phase. The Matchete v0.6.0 FeynmanRules return shape
has not been tested on the user's install; the output will determine how to
assemble the signed amplitude without guessing the spinor-flow convention.

## PowerShell

```
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$dir = 'output\general_observable\lo_fermion_boundary'
Expand-Archive -LiteralPath "$HOME\Downloads\T3_Native_Feynman_Rules_A_m4.zip" -DestinationPath . -Force
python -X utf8 "$dir\test_native_feynman_rules.py"
python -X utf8 "$dir\native_feynman_rules_check.py"
$env:T3_EXPERIMENTAL_GENERAL_FERMIONS = '1'
$env:T3_EXPERIMENTAL_GENERAL_QUARTICS = '1'
$env:T3_EXPERIMENTAL_GENERAL_QUADRATICS = '1'
$env:T3_EXPERIMENTAL_FLAVOR_MASSES = '1'
& wolframscript -file "$dir\native_feynman_rules_audit.wl" *> "$dir\normalization_probes\native_feynman_a_m4_run.log"
Get-Content "$dir\normalization_probes\native_feynman_a_m4_run.log" -Tail 15
$report = Get-ChildItem "$dir\normalization_probes\native_feynman_a_m4_*.json" | Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName
python -X utf8 "$dir\native_feynman_rules_check.py" --report $report
```

If blocked, the JSON and log identify the stage; do not label the boundary as
physically normalized. Share the generated `native_feynman_a_m4_*.json`
(or ZIP it if too large). No other model is changed.
