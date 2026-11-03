# A4 — Exact PowerShell commands (Windows, no overwrite)

All paths are **project-relative** under `C:\thesistemp\T3neutrinomassmodelpipeline`; the ZIP itself can be placed anywhere. The ZIP contains **only additive files** under `Numerical/lo_certification/physical_normalisation/parallel_A_spinor/`. Nothing writes to GitHub or modifies the production Lagrangian, matching, RGE, or provisional one-loop engine.

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = 'C:\path\where\you\saved\TrackA_B1_CommonSpinor.zip'
$track = 'Numerical\lo_certification\physical_normalisation\parallel_A_spinor'
if (Test-Path (Join-Path $repo $track)) { throw 'Refusing to overwrite an existing Track A directory' }
Expand-Archive -LiteralPath $zip -DestinationPath $repo
Set-Location $repo
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss_ffff'
$out = "output\general_observable\physical_normalisation\parallel_A_spinor_$stamp"
if (Test-Path $out) { throw 'Refusing to overwrite output directory' }
New-Item -ItemType Directory -Path $out -Force:$false | Out-Null

# A1: exact exterior-algebra spinor witness (no production imports).
python "$track\weyl_algebra_check.py" --output "$out\exact_weyl_algebra.json"
if ($LASTEXITCODE -ne 0) { throw 'A1 exact Weyl check failed' }

# A2: NEW installed-Matchete rule-level check. Needs Matchete with FeynmanRules,
# typically v0.6+; it does not repeat scalar/CG/mass witness batches.
$rules = "$out\native_feynman_rules.json"
wolframscript -file "$track\native_rule_probe.wl" $rules
if ($LASTEXITCODE -ne 0) { throw 'A2 native FeynmanRules unavailable or native rule/negative control blocked' }

# A3: source-native witness-gated complex UV calculations.
# Requires the three *already produced* native scalar/CG JSONs and the existing
# B1 spinor action JSON at the exact project paths listed in complex_three_models.py.
python "$track\complex_three_models.py" `
    --mode native `
    --repo-root . `
    --rule-witness $rules `
    --output-dir "$out\complex_native"
if ($LASTEXITCODE -ne 0) { throw 'A3 native-source provenance or numerical checks failed' }

# Print the actual statuses; do not upgrade physical certification.
Get-Content "$out\complex_native\three_model_complex_results.json" -Raw |
    ConvertFrom-Json |
    Select-Object mode, native_GammaCC_relative_phase_fixed, physical_source_certified, certified_count_of_16
```

`A2` may **block** even though the earlier installed-Matchete mass witnesses succeeded: public `FeynmanRules` arrived in Matchete v0.6.0 and the installed version/capability has not been examined in this environment. A successful `FeynmanRules` invocation still does **not** establish a relative `GammaCC`↔Weyl component phase; the script records the exact Feynman rule text for manual spinor-basis analysis. The missing numerical `GammaCC` basis is a **physics obstruction**, not merely a package-installation error.

For an independent **non-native** replay of the included mathematical calculations, with no installed Matchete and no native-file prerequisites:

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss_ffff'
$offline = "output\general_observable\physical_normalisation\parallel_A_spinor_offline_$stamp"
if (Test-Path $offline) { throw 'Output already exists' }
python 'Numerical\lo_certification\physical_normalisation\parallel_A_spinor\complex_three_models.py' `
    --mode audited-constants --repo-root . --output-dir $offline
if ($LASTEXITCODE -ne 0) { throw 'Independent declared-convention numerical check failed' }
```

Dependencies for the Python executables: `python`, `numpy`, `mpmath`, `sympy`. For native A2: `wolframscript`, `Matchete` loaded via `Needs["Matchete`"]`, and a version exposing `FeynmanRules`. No source changes are needed.

The A2 probe exports rich but **non-overwriting** native `InputForm` records. The output can be large. Do not overwrite a prior run; change the timestamp. If any native output is absent or failed, label `A2 BLOCKED`, `A3 physical source NOT CERTIFIED`, and keep `0/16`.
