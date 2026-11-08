# Workstream B — B/C Higgs-mass RGBeta WolframScript fix v2

**Scope:** replace `RGE/running/rgbeta/RunT3BCFullScalarMassRGBeta.wl` only; add a source-level regression test. No model or production RGE code is changed. Both B(-1) and C(-1), UV/EFT1 remain covered.

## What happened

The first live `B/UV` run printed `FrontEndObject::notavail` and then failed the `A required one-loop mass beta is unresolved` gate, so none of the four JSON exports existed. The Python verifier was correct to reject missing evidence.

The v1 runner used `Check[BetaTerm[sym,1],$Failed]`, which returns `$Failed` after *any* message generated during evaluation, even if the beta expression itself is valid. The source does not establish whether the front-end message actually caused that rejection; this is the **most likely mechanism**. `UsingFrontEnd` has already worked for the project's native Matchete probes on the user's machine.

## Changes in v2

1. Evaluate the complete experimental runner within `UsingFrontEnd`.
2. Replace each broad message-sensitive `Check[BetaTerm[...]]` with `CheckAbort[Quiet[BetaTerm[...]],$Aborted]`. This does not convert a genuine `$Failed`, unevaluated `BetaTerm`, `Missing`, `Indeterminate`, or `ComplexInfinity` to success.
3. Print the coupling symbol before requesting each mass beta. If any remain unresolved, save a `Blocked` JSON report containing `beta_diagnostics` for each symbol, then exit nonzero. The successful JSON format and the existing `VerifyT3BCFullScalarMassClosure.py` contract are unchanged.

**No live beta has been certified by local Python tests.** The standalone v2 runner still needs to be executed using the installed RGBeta/Wolfram stack.

## Run in PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath '.\Workstream_B_BC_Higgs_Mass_FrontEnd_Fix_v2.zip' -DestinationPath . -Force
python -X utf8 -m pytest -q tests\test_t3_bc_full_scalar_mass_closure.py tests\test_t3_bc_full_scalar_mass_runner_v2.py

$root = 'output\audits\workstream_b_bc_full_mass'
foreach ($cls in @('B', 'C')) {
    foreach ($stage in @('UV', 'EFT1')) {
        $name = "${cls}_m1_$($stage.ToLower())_full_mass_v1.json"
        wolframscript -file RGE\running\rgbeta\RunT3BCFullScalarMassRGBeta.wl $stage $cls "$root\$name"
        if ($LASTEXITCODE -ne 0) {
            Write-Host "Blocked: $cls/$stage (exit $LASTEXITCODE)"
            if (Test-Path "$root\$name") {
                python -c "import json,sys; x=json.load(open(sys.argv[1])); print(json.dumps({'status':x.get('status'),'error':x.get('error'),'beta_diagnostics':x.get('beta_diagnostics')},indent=2))" "$root\$name"
            }
            throw 'Wolfram beta export blocked'
        }
    }
}
python -X utf8 -m scripts.VerifyT3BCFullScalarMassClosure --directory $root --output "$root\closure_v1.json"
```

Note: reports will not be overwritten. The failed original `B/UV` v1 attempt wrote no file according to the user output. If a report already exists, choose a fresh **directory** and keep the original basename pattern ending in `_full_mass_v1.json` because the existing verifier expects it.

The closure verifier checks the portal coefficients and benchmark mass-only evolution. It does **not** certify the Matchete/RGBeta field phase, threshold matching, or full running of dimensionless couplings.
