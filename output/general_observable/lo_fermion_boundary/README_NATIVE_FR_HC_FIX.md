# A_-4 native Matchete Feynman-rule input fix

Replaces the complete `native_feynman_rules_audit.wl` file in
`output/general_observable/lo_fermion_boundary/` only.

## Verified bug

The v1 diagnostic passed only `contact`, the single holomorphic `Proj[-1]`
LLSS source term, to `FeynmanRules`. Matchete v0.6.0 rejected this input
with `CheckLagrangian::Hermiticity` and reported an empty EFT1 contact vertex.

The fixed diagnostic passes

```mathematica
contactWithHC = RelabelIndices[PlusHc[contact]];
eftRules = CheckAbort[UsingFrontEnd[FeynmanRules[contactWithHC]], $Aborted];
```

The full term, including complex Yukawa coefficients, is conjugated by
Matchete's native `PlusHc`, not reconstructed by hand. Both `Proj[-1]` and
`Proj[1]` sectors are required in the resulting expression. The
holomorphic Wilson coefficient **is not doubled**; the HC is present only
to supply a valid input Lagrangian to the Feynman-rule routine.

If extraction still fails, the report now captures the original PL contact,
Hermitian candidate, UV vertices, and attempted EFT rule in the new JSON.

## PowerShell

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$dir = 'output\general_observable\lo_fermion_boundary'
Expand-Archive -LiteralPath "$HOME\Downloads\T3_Native_FR_Hermiticity_Fix.zip" -DestinationPath . -Force
python -X utf8 "$dir\test_native_feynman_rules.py"
& wolframscript -file "$dir\native_feynman_rules_audit.wl" *> "$dir\normalization_probes\native_feynman_a_m4_run_hc.log"
Get-Content "$dir\normalization_probes\native_feynman_a_m4_run_hc.log" -Tail 20
$report = Get-ChildItem "$dir\normalization_probes\native_feynman_a_m4_*.json" |
  Sort-Object LastWriteTime -Descending |
  Select-Object -First 1 -ExpandProperty FullName
python -X utf8 "$dir\native_feynman_rules_check.py" --report $report
```

Requires the four `T3_EXPERIMENTAL_*` variables already set to `1`.
No source or certified numerical boundary is rewritten. A successful vertex
extraction is NOT yet a connected UV amplitude or an absolute spinor
normalization certificate.
