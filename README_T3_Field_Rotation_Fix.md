# T3 self-bilinear: Matchete `Field[...]` rotation fix

## Exact project destination

- `Lagrangian/RunT3GeneralSelfMassBasis.wl` — **complete replacement file**.

This revision fixes the exit-code-6 rotation guard observed in A(-2), where
Matchete's fully evaluated fields are of the form
`Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]`, not
`Bar[NewScalar2[d$$1]]`. It rotates both `Field[...]` and `Bar[Field[...]]` in
the complete general interactions, preserving all Yukawas and quartics. The
holomorphic bilinear remains encoded by the separate split masses. The rest
of the runner and the existing `RunT3GeneralMatcheteMatching.py` are unchanged.

The runner checks that no original field or holomorphic mass insertion remains,
and that both mass eigenfields are present. The old diagnostic JSON is removed
at the start of each new attempt to prevent stale failures being mistaken for
new ones.

## Install (PowerShell, project root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive "$HOME\Downloads\T3_Self_Bilinear_Field_Rotation_Fix.zip" -DestinationPath . -Force
```

## Run one symbolic matching test

```powershell
python -X utf8 -m scripts.RunT3GeneralMatcheteMatching --model A_-2
```

## Inspect output

```powershell
$dir = "output\general_matchete_matching\A_m2"
Get-Content "$dir\general_matching_status.json" -Raw
Get-Content "$dir\general_matching.log" -Tail 100
if (Test-Path "$dir\self_mass_basis_rotation_debug.json") {
    Get-Content "$dir\self_mass_basis_rotation_debug.json" -Raw
}
```

**Expected rotation-stage diagnostics** if the field rules work:

```text
ROTATION FIELD WRAPPERS: before=<positive count>; after=0
ROTATION CHECK: old field before RelabelIndices=0; after=0; extra mass before=0; after=0
```

This does **not** predict that Matchete will complete. Exit code 7 denotes a
rotated `CheckLagrangian` failure; code 8 denotes matching failure; code 9
denotes incomplete/missing Weinberg coefficient. Those will require separate
fixes based on the live log.

## Verification / remaining limitations

- Verified with mock `Field[...]/Bar[Field[...]]` Wolfram Language expressions:
  the original field vanishes and both real eigenfields are present for S1
  and S2 rotations.
- The complete source passed offline bracket/comment/string balance checks.
- Not verified using the user's local Matchete installation.
- Real-positive bilinear convention only; complex phases, physical
  Matchete/RGBeta normalization, and sequential EFT Wilson matching are not
  implemented/certified by this fix.
- **Do not use new fitted chi-square values until actual general Matchete
  mass-basis matching succeeds.**
