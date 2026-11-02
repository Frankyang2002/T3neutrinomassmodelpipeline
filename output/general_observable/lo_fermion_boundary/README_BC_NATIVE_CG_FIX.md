# B_-1/C_-1 native CG extractor: frontend and kernel-state fix

This replaces only the native B/C CG extractor, its Python validator and
its offline tests. All files reside under
`output/general_observable/lo_fermion_boundary/`. There is **no** change to
the original 16-model boundary or the Matchete model/matching source.

## Reason for change

The previous extractor called `UsingFrontEnd[Needs["Matchete`"]]` but ran
`ResetAll[]`, `LoadModel["SM"]`, and `PhysicalSU2InvariantBasis` **after** leaving
the frontend scope. The supplied log shows `FrontEndObject::notavail` and
a subsequent failure during C_-1 initialization. The extractor now puts
all Matchete initialization and native-invariant work inside `UsingFrontEnd`
and executes **one model per fresh wolframscript kernel**. This separates
C_-1 from any B_-1 registry/cache state. A genuine frontend failure remains
fail-closed with a per-model `Blocked` JSON report.

The previous Python `--report $report` error was secondary: no successful
combined JSON existed so PowerShell passed an empty argument.

## Installed paths

- `output/general_observable/lo_fermion_boundary/native_bc_yukawa_cg_probe.wl`
- `output/general_observable/lo_fermion_boundary/check_bc_native_cg_products.py`
- `output/general_observable/lo_fermion_boundary/test_bc_native_cg_products.py`
- `output/general_observable/lo_fermion_boundary/README_BC_NATIVE_CG_FIX.md`

## Run (PowerShell, from project root)

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$dir = 'output\general_observable\lo_fermion_boundary'
$zip = "$HOME\Downloads\T3_BC_Native_CG_FrontEnd_Fix.zip"
Test-Path -LiteralPath $zip
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

python -X utf8 "$dir\test_bc_native_cg_products.py"
if ($LASTEXITCODE -ne 0) { throw 'Python tests failed' }

foreach ($model in @('B_-1','C_-1')) {
    $tag = $model.Replace('-', 'm').Replace('+', 'p')
    $log = "$dir\normalization_probes\bc_native_yukawa_cg_${tag}_run.log"
    & wolframscript -file "$dir\native_bc_yukawa_cg_probe.wl" $model *> $log
    Get-Content -LiteralPath $log -Tail 12
    if ($LASTEXITCODE -ne 0) { throw "$model native Matchete extractor failed; inspect $log" }
}

$folder = "$dir\normalization_probes"
$b = Get-ChildItem $folder -Filter 'bc_native_yukawa_cg_B_m1_*.json' |
     Sort-Object LastWriteTime -Descending |
     Select-Object -First 1 -ExpandProperty FullName
$c = Get-ChildItem $folder -Filter 'bc_native_yukawa_cg_C_m1_*.json' |
     Sort-Object LastWriteTime -Descending |
     Select-Object -First 1 -ExpandProperty FullName
if (-not $b -or -not $c) { throw 'Missing B/C native CG report(s)' }

python -X utf8 "$dir\check_bc_native_cg_products.py" --report-b "$b" --report-c "$c" --check-only
if ($LASTEXITCODE -ne 0) { throw 'Native CG contraction check failed' }
python -X utf8 "$dir\check_bc_native_cg_products.py" --report-b "$b" --report-c "$c"
if ($LASTEXITCODE -ne 0) { throw 'Native CG report validation failed' }
```

The Python audit reads the existing pinned B/C source, Matchete-generated
matched CG registry, and formal release. It checks Yukawa tensor orientations,
B_-1 identity matrices, and the C_-1 predicted decomposition

`K[i,j,a,b] = Q[i,j,a,b]/sqrt(3) + sqrt(3)/2 * eps[i,j] * Bridge[b,a]`

for the 11/12/22 Yukawa pairs. The native `-alpha` y2 sign is separately
verified. This is a **test target**, not a forced normalization; disagreement
fails with exit code 2. A pass does *not* independently evaluate the full
Majorana Feynman-rule amplitude or certify the absolute `GammaCC`/master-Weyl
phase. The sixteenth-model release remains unchanged.

## Status

The 11 Python tests have run successfully in the preparation environment.
The corrected Mathematica extractor has **not** been executed there because
Wolfram/Matchete is only available on the user's machine. The displayed
frontend failure from the older script is not itself evidence of a CG mismatch.
