# A(-4) general Higgs–triplet portal: production source correction

## Verified evidence

- A(-4) full-general UV Lagrangian builds with both `HiggsPortal2Inv1` and `HiggsPortal2Inv2`.
- The original tree `RawEFT` contains both; `GreensSimplify` removes `lambdaH2Inv1` (but retains `lambdaH2Inv2`).
- The user's native CG test established `Inv1` is **anti-Hermitian**, `Inv2` is **Hermitian**, and the independent CG rank is 2.
- The user's in-memory test succeeded with `I lambdaH2Inv1` **inside** `PlusHc`, and both independent portals survived the complete F-only Matchete tree matching.

## Production change

Full replacement `Lagrangian/interactions/ScalarPotential.wl` retains the existing general potential generator and other operator families, but assigns the Hermitian phase **from the native CG tensor** to each Higgs–scalar portal:

```
T^HC_{ijab} = Conjugate[T_{ji ba}]
T == T^HC    -> (1/2) PlusHc[lambda O]
T == -T^HC   -> (1/2) PlusHc[I lambda O]
```

For singlet scalars, only the Higgs pair is swapped. The native tensor must be nonzero with the expected dimensions. The generator fails closed if its tensor is neither Hermitian nor anti-Hermitian under this map. It does *not* assume that the index of the anti-Hermitian invariant is universally `Inv1` across all models.

For A(-4), the change reproduces the validated trial mechanism:

```
1/2 PlusHc[I lambdaH2Inv1[] Hdag H S2dag S2 CG[Inv1]]
```

The explicit `I` is required because the original tensor is anti-Hermitian. It is **not** a certified numerical sign/normalisation relation to RGBeta's `lambdaH2Adj`.

## Included files (paths relative to repository root)

- `Lagrangian/interactions/ScalarPotential.wl` — **FULL replacement** (the only production file changed).
- `Lagrangian/TestAminus4HermitianPortalProduction.wl` — additive fresh-kernel regression; uses the installed builder with **no override**; checks both portals after `RunT3Matching`.
- `EFT1_A_M4_GENERAL_PORTAL_PRODUCTION_FIX.md` — this guide.

The source base is the latest GitHub `Lagrangian/interactions/ScalarPotential.wl` inspected on 2026-11-09, SHA `58f5bbef2688f9e87a1b8992bca11076ee254437`. No GitHub files were modified. Only the A(-4) Matchete trial is physically exercised; this does **not** certify all 16 models.

## Install from PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Copy-Item .\Lagrangian\interactions\ScalarPotential.wl .\Lagrangian\interactions\ScalarPotential.wl.before_portal_fix
Expand-Archive -LiteralPath "$HOME\Downloads\T3_EFT1_Aminus4_Hermitian_Portal_Production.zip" -DestinationPath . -Force
$env:T3_EXPERIMENTAL_GENERAL_FERMIONS = '1'
$env:T3_EXPERIMENTAL_GENERAL_QUARTICS = '1'
$env:T3_EXPERIMENTAL_GENERAL_QUADRATICS = '1'
$env:T3_EXPERIMENTAL_FLAVOR_MASSES = '1'
wolframscript -file Lagrangian\TestAminus4HermitianPortalProduction.wl `
  output\eft1_general_lagrangian_trace\A_m4_production_regression_01
Get-Content 'output\eft1_general_lagrangian_trace\A_m4_production_regression_01\production_portal_regression.json' -Raw
```

Expected success status: `BothGeneralAminus4PortalsRetained_InventoryOnly`.

## Regenerate the complete matched A(-4) F-first source **to a fresh folder**

Only after the regression succeeds:

```powershell
$out = 'output\general_observable\fermion_first_hermitian_portal_v1\A_m4'
New-Item -ItemType Directory -Path $out -Force | Out-Null
& wolframscript -file Lagrangian\RunT3GeneralFThreshold.wl $out A_-4 0 *>&1 |
  Tee-Object -FilePath (Join-Path $out 'native.log')
if ($LASTEXITCODE -ne 0) { throw 'F-first Matchete regeneration failed' }
python -X utf8 -m scripts.AuditAminus4GeneralLagrangian `
  --source-dir $out `
  --output 'output\eft1_general_lagrangian_audit\A_m4_hermitian_v1.json'
```

Expected inventory status after successful regeneration: `GeneralRenormalizableEFT1InventoryPresent`, and no missing quartics. That status only certifies the **matched renormalisable operator inventory**, not the physical coupling dictionary.

**Do not rerun the Stage B saved-benchmark ODE expecting the corrected portal automatically.** Stage B still reads the *old* `batch16_v1` source, its matching hash, and saved renormalisable boundary. A fresh, provenance-consistent benchmark must be constructed before physically combining the regenerated source and the Wilson running. Never mix F-matching and beta/coupling files from different benchmark points.

Remaining: exact Matchete ↔ RGBeta scalar CG and coupling normalisation, `GammaCC` ↔ Weyl spinor factor, full dimension-five closure, simultaneous **both-scalar** matching at M_S, SMEFT handoff. None is certified by this correction.
