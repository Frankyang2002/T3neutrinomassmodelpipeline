# Workstream C: Matchete v0.6 NCM-versus-** operator identity (v11.2)

## Verified cause of the v11.1 failure

The user's v11.1 local report reads:

- `direct_equivalence_status = "NotProvenByNativeSimplifier"`.
- Installed Matchete and built-in SM both loaded successfully.
- The reference spinor is `NCM[Transp[l], DiracProduct[GammaCC,Proj[-1]]] ** NCM[DiracProduct[Proj[-1]], l]`.
- The native spinor is `NCM[Transp[l], DiracProduct[GammaCC,Proj[-1]], l]`.

These are **different symbolic heads**. Matchete's official v0.4.0 changelog explicitly states it replaced `NonCommutativeMultiply` (`**`) with `NCM` for spinor products, retaining `\[CenterDot]` as the corresponding infix operator. The installed `Models/SMEFT_Warsaw.m` source, supplied via `GammaCC_source_evidence.zip`, still contains an old `**` expression for `Lagd5`. Thus the v11.1 comparison was not a valid test of the desired spinor equality; its failure does **not** establish a sign difference.

## Changes

One full replacement file:

`Lagrangian/ProbeT3BCWarsawDirectIdentity.experimental.wl`

New accompanying explanation:

`docs/WorkstreamC_WarsawDirectIdentity_v11_2_README.md`

v11.2 preserves the v11.1 frontend/SM initialization, the original operator ordering, and the Python v11 checker. The only physics-relevant change is to construct the installed Warsaw spinor product using `NCM[Bar@CConj@l[j,p], l[m,r]]`, instead of its legacy `Bar@CConj@l[j,p] ** l[m,r]` spelling. The probe exports both spellings for verification and explicitly tests both positive and negative operator equality. It does not substitute a sign by hand.

The old `**` construction remains present solely as *diagnostic evidence*. It is **never** used as the reference for the equality test.

## Run (PowerShell at repository root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "$HOME\Downloads\T3_WorkstreamC_WarsawDirectIdentity_v11_2.zip" `
  -DestinationPath . -Force

$dir = "output\audits\t3_bc_y_warsaw_direct_identity_v11"
wolframscript -file Lagrangian\ProbeT3BCWarsawDirectIdentity.experimental.wl `
  "$dir\native_identity.json"
$probeExit = $LASTEXITCODE
Get-Content "$dir\native_identity.json"

if ($probeExit -eq 0) {
  python scripts\ResolveT3BCWeinbergOperator.experimental.py `
    --matchete-identity "$dir\native_identity.json" `
    --output "$dir\summary_live.json"
}
```

## Interpretation

- `ExactPlusEquality`: the Matchete NCM reference equals the native PL operator with **positive relative sign** under the installed conventions. The original v11 Python check can then confirm the chiral orientation and the factor converting symmetric Warsaw coefficients to conventional `1/2`-normalized `C5`.
- `ExactMinusEquality`: opposite sign; **fail**, do not alter source conventions automatically.
- `NotProvenByNativeSimplifier`: no identity established; the report includes both reduced `InputForm` strings and their difference.
- `UnexpectedNoncommutativeHead` or `DegenerateZeroOperator`: **fail**; source or algebra assumptions changed.

The explicit use of `NCM` follows the current Matchete algebra, but no claim of successful **local** Wolfram evaluation is made until this script is run in the user's installation.

This operator check is **not** a proof of full EFT1/EFT2 equivalence, EFT counterterm subtraction, physical renormalized Y/X finite matching, or physical `C5`.

Official Matchete changelog: <https://matchete.gitlab.io/changelog.html> (v0.4.0 entry).
