# Direct UV one-loop: native scalar-portal normalization audit

## What is being checked

The source `Lagrangian/interactions/ScalarPotential.wl` constructs independent
`HiggsPortal{1,2}Inv{k}` invariants and native Matchete couplings
`lambdaH{1,2}Inv{k}`. The source term is `1/2 PlusHc[lambda * Hdag H * Sdag S * CG]`.
For `H^0=v/sqrt(2)`, its contribution to the quadratic scalar potential is

```
Delta A^(i)_{st} = -v^2/2 * Sum_k lambdaHiInvk * P^(i,k)_{st}
P^(i,k) = HermitianPart[T^(i,k)_{Bar H0,H0,Bar S,S}].
```

The frozen numerical references contain `lambdaH1`, `lambdaH2`, which are **not**
a demonstrated physical dictionary for these distinct native Matchete couplings.
The existing provisional adapter instead assumes `Delta A = -v^2/2 * lambdaHi * I`.
The audit exposes the Matchete tensor projectors, and reports any disagreement
under the **explicitly conditional** hypothesis `lambdaHiInv1=lambdaHi`, with
all other native portal invariants zero. It does not establish that hypothesis.

## Files added or replaced

- `Lagrangian/ExportT3DirectUVHiggsPortalCG.wl` (new): export native registered
  portal tensors for one model.
- `scripts/AuditT3DirectUVPortalBasis.py` (new): project all 16 into the
  charge basis verified earlier, with source-contract guards and provenance.
- `scripts/RunT3DirectOneLoop.py` (full replacement): fix the **misleading**
  terminal summary, previously printing "Source-certified ... 16/16" for
  16 *uncertified* provisional evaluations. Numerical algorithms are unchanged.
- `tests/python/test_t3_direct_one_loop_native_portal.py` (new): pure guards.
- `tests/python/test_t3_direct_one_loop_runner_label.py` (new): label regression.

## Location for newly generated reports

All generated files stay below
`output/general_observable/lo_direct_one_loop/`:

- `source_portal_cg_all16/T3_.../native_higgs_portal_registry.json`
- `source_portal_cg_all16/T3_.../export.log`
- `source_audits/native_higgs_portal_basis_audit.json`

The existing `provisional/` and `provisional_inputs/` are **never changed** by
the new exporter/audit. No EFT code or GitHub repository is modified.

## PowerShell (repository root)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
python -m unittest discover -s tests/python -p "test_t3_direct_one_loop*.py" -v

$root = (Get-Location).Path
$base = Join-Path $root "output\general_observable\lo_direct_one_loop"
$portalRoot = Join-Path $base "source_portal_cg_all16"
$auditDir = Join-Path $base "source_audits"
New-Item -ItemType Directory -Force -Path $portalRoot, $auditDir | Out-Null
$configs = @(Get-ChildItem "Numerical\benchmarks\fitted_reference\T3_*.json")
if ($configs.Count -ne 16) { throw "Need exactly 16 frozen benchmarks" }
foreach ($file in $configs) {
    $cfg = Get-Content $file.FullName -Raw | ConvertFrom-Json
    $r = $cfg.representation
    $modelDir = Join-Path $portalRoot $file.BaseName
    New-Item -ItemType Directory -Force -Path $modelDir | Out-Null
    $dest = Join-Path $modelDir "native_higgs_portal_registry.json"
    $log = Join-Path $modelDir "export.log"
    & wolframscript -file "Lagrangian/ExportT3DirectUVHiggsPortalCG.wl" `
        "$([int]$r.d_s1)" "$([int]$r.d_s2)" "$([int]$r.d_f)" `
        "($([int]$r.alpha))" "$dest" *> "$log"
    if ($LASTEXITCODE -ne 0) {
        Get-Content $log -Tail 60
        throw "Native portal CG export failed: $($file.BaseName)"
    }
    $export = Get-Content $dest -Raw | ConvertFrom-Json
    if ($export.status -ne "Success") { throw "Incomplete portal: $($file.BaseName)" }
    Write-Host "Native portal: $($file.BaseName) PASS"
}

python scripts/AuditT3DirectUVPortalBasis.py `
    --project-root . `
    --portal-root "$portalRoot" `
    --audit-dir "$auditDir" `
    --output (Join-Path $auditDir "native_higgs_portal_basis_audit.json")
if ($LASTEXITCODE -ne 0) { throw "Native portal basis audit failed" }
```

## What success does and does not prove

If all 16 export and the audit succeeds, the *native* Higgs-portal tensors and
their neutral-Higgs, charge-diagonal components are known, in the tested
candidate charge basis. The output does **not** certify that the frozen
`lambdaH1,lambdaH2` parameters equal Matchete's `lambdaH1Invk,lambdaH2Invk`,
nor that the provisional neutrino mass matrices are physically normalized.
Spinor charge conjugation, fermion phase transport, optional Yukawa couplings,
and any extra scalar quadratic interactions remain distinct workstreams.

Do not reinterpret a conditional difference as a proven bug in the frozen
RGBeta parameterization. It is a gap in the presently unverified dictionary.
