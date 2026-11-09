# T3 A(-4) Stage B — full-general Lagrangian source gate

Source review: public `Frankyang2002/T3neutrinomassmodelpipeline` default-branch commit
`ea93a730812cc6a288ae7eba363874a2166ea6c6` (read-only), 2026-11-09.

## Verified distinction

- The saved F-only Matchete producer `Lagrangian/RunT3GeneralFThreshold.wl`
  explicitly requires the four `T3_EXPERIMENTAL_GENERAL_*` flags and calls
  `BuildT3Lagrangian` with the full interaction inventory. The saved native log
  records both `HiggsPortal2Inv1` and `HiggsPortal2Inv2` as individually valid
  and reports **two independent** H-S2 portal invariants.
- The **matched** `eft1_after_F_tree_inputform.txt` does not contain
  `Coupling[lambdaH2Inv1,...]`; it contains only `lambdaH2Inv2` of the two
  independently validated portals. Other expected potential couplings are
  present. This is an unresolved loss/redefinition in the saved matching
  output, **not** proof that the UV Lagrangian lacked the portal.
- A(-4) needs no additional *new* general quartics beyond the established
  families: `Lagrangian/interactions/GeneralScalarVerifiedQuartics.wl`
  lists `1_3_-4 -> {}`. That does **not** remove either of the two ordinary
  H-S2 portal invariants.
- The archived EFT1 RGBeta payload has 17 beta functions including both
  `lambdaH2` and `lambdaH2Adj`. The master Python reference quartic basis
  includes both, but the Matchete->RGBeta scalar tensor and coupling map has
  **not** been certified. The existing `lambdaH2Adj` source-normalization
  comparison (independent real-tensor reference `1` versus RGBeta `3`)
  remains unresolved.
- `Numerical/running/T3EFT1WilsonODETransport.py` uses the full 7-term
  `psi^2 phi^2` master RGE but currently an **assumed identity** scalar map
  and an unresolved `GammaCC`-to-Weyl overall factor. Other possible
  dimension-five source classes have not been proven absent.

**Conclusion:** prior Stage B A(-4) C5 is a useful *formal* numerical ODE test,
not a source-complete general-UV->EFT1 physical prediction. The strict new
runner intentionally stops at the source gate on this saved data. It does not
invalidate the earlier numerical convergence test, but that test is only
within an unverified input convention.

## Files included in this ZIP

1. `scripts/AuditAminus4GeneralLagrangian.py` — new read-only strict audit.
2. `scripts/RunAminus4SavedEFT1Pilot.py` — **complete replacement** of the
   previous Stage B runner. Adds the source-inventory check *before* creating
   the run directory. Saves the audit in `input_validation.json` if it passes.
3. `tests/test_a_minus4_general_lagrangian_audit.py` — four synthetic unit
   tests: complete inventory, missing portal, missing RGBeta coupling and
   false portal-invariant claim.

No existing source outputs or saved numerical results are changed. The user
must regenerate/resolve the Matchete matching source before the strict runner
will permit further runs.

## PowerShell

From the repository root:

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = Join-Path $HOME 'Downloads\T3_EFT1_General_Lagrangian_Gate.zip'
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p 'test_a_minus4_general_lagrangian_audit.py' -v
python -X utf8 -m unittest discover -s tests -p 'test_a_minus4_saved_eft1_pilot.py' -v
python -X utf8 -m scripts.AuditAminus4GeneralLagrangian
```

**Expected with the current saved inputs:** audit status
`Blocked_GeneralEFT1InventoryNotEstablished`, missing
`lambdaH2Inv1`, process exit code **2**. PowerShell will display an error
exit status; this is deliberate. The audit still writes:

`output/eft1_general_lagrangian_audit/A_m4.json`

Inspect it:

```powershell
Get-Content 'output\eft1_general_lagrangian_audit\A_m4.json' -Raw
```

An attempt to rerun `RunAminus4SavedEFT1Pilot` with the current source now
raises a `SavedPilotError` before creating an output directory; previous
`initial_ode` and `convergence_15` outputs remain untouched.

## Required source-level correction / next numerical step

Trace `lambdaH2Inv1` through these steps **in the same fresh Wolfram kernel**:

1. `Lagrangian/model/LagrangianBuilder.wl`: ensure `LUV` contains *both*
   `lambdaH2Inv1` and `lambdaH2Inv2` and that their CG tensors are linearly
   independent after SU(2) index contraction.
2. `Lagrangian/RunMatching.wl`: identify the exact stage where one independent
   invariant disappears (tree `Match`, `GreensSimplify`, `EOMSimplify`,
   `ReplaceEffectiveCouplings`, or export). Do not assert its cause before
   this diagnostic.
3. Re-export the corrected general EFT1 matched source; establish the explicit
   Matchete->RGBeta quartic/basis map. Verify all 9 expected source quartic
   **coupling families** and all 17 RGBeta EFT1 beta keys.
4. Rebuild the saved fermion-first boundary and same-run EFT1 trajectory before
   rerunning the Stage B Wilson ODE. The old source hash will necessarily
   change; do not patch it or reuse its C5.

This gate checks operator *inventory* only. Even a passing audit does not
certify GammaCC/Weyl normalization, physical CG basis, dimension-five
operator closure, finite scalar matching, or SMEFT evolution.
