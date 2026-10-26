# Stage 22 — Local B/C reused quartic SU(2) correction

## Source and scope

Full replacement `RGE/running/rgbeta/T3RGBetaModel.wl` is reconstructed from the
current GitHub `main` blob `584ccb3ada2a20b980442533bf0dac5160d0ff10`
(confirmed byte-for-byte with `git hash-object`). The **only change** is

```wolfram
(* Before *)
eps[SU2L @ fund, #1, #3] * eps[SU2L @ fund, #2, #4] &
(* After *)
del[SU2L @ fund, #1, #3] * eps[SU2L @ fund, #2, #4] &
```

inside `lambdaHHdagS1barS2barCross`, field order
`{H, Bar @ H, Bar @ S1, Bar @ S2}`. Other RGBeta model interactions unchanged.
The original live RGBeta native covariance export registered two SU(2) structures
in B_-1 and C_-1. It found generator nonzero residual counts [0,0,0]
(new `lambdaT3GeneralQ110101I1`) and [4,0,4] (reused coupling).
The proposed delta-epsilon tensor passes all three generators in exact Python
analysis, and the two RGBeta geometry columns have exact basis matrix

    R = M A,  A = [[1/2, 1], [-sqrt(3)/2, 0]], det A = sqrt(3)/2.

This **does not** prove RGBeta component-basis identification outside this
block, establish signed scalar-potential prefactors, or certify the physical
Matchete–RGBeta coupling map. Stage 22 continues to FAIL CLOSED.

## Changed/generated full files

- `RGE/running/rgbeta/T3RGBetaModel.wl`: full production replacement; one contraction changed.
- `scripts/Stage22RGBeta110101Covariance.py`: supports old and corrected native exports, never silently accepts other tensors.
- `scripts/Stage22RGBetaAlignedGeometry.py`: supports historical obstruction or exact corrected 2D block.
- `scripts/VerifyGeneralT3Stage22Live110101Fix.py`: strict check of independently generated live Wolfram native tensor components and covariance.
- `tests/python/test_general_t3_stage22_110101_covariance.py`: pre/post and synthetic correction tests.
- `tests/python/test_general_t3_stage22_aligned_geometry.py`: old and synthetic corrected block tests.
- `tests/python/test_general_t3_stage22_explicit_basis.py`: expected pre/post inventory totals.
- `tests/python/test_general_t3_stage22_live_110101_fix.py`: live-Wolfram verifier tests.
- `docs/GENERAL_T3_STAGE22_110101_LOCAL_CORRECTION.md`: this document.

## Windows PowerShell

From repository root, back up any current local model file before extracting:

```powershell
Copy-Item .\RGE\running\rgbeta\T3RGBetaModel.wl .\RGE\running\rgbeta\T3RGBetaModel.wl.pre-stage22.bak
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage22_110101_LocalFix.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_stage22.py tests\python\test_general_t3_stage22_explicit_basis.py tests\python\test_general_t3_stage22_aligned_geometry.py tests\python\test_general_t3_stage22_110101_covariance.py
```

Then run the existing Wolfram exporter to ensure RGBeta registered the
corrected function in **both** models:

```powershell
wolframscript -file .\scripts\ExportGeneralT3Stage22SU2NativeCovariance.wl output/general_t3_stage22/stage22_110101_wolfram_native_covariance.json
python -m scripts.VerifyGeneralT3Stage22Live110101Fix
python -m pytest -q tests\python\test_general_t3_stage22_live_110101_fix.py
```

After successful live audit, regenerate Stage 22 declarations and invariants
without repeating the Matchete export:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage22.ps1 `
    -ExportRGBetaDeclarations -NativeDiagnosticsOnly
python -m scripts.Stage22RGBeta110101Covariance
python -m scripts.Stage22RGBetaAlignedGeometry
python -m scripts.Stage22RGBetaExplicitBasis
```

Expected post-fix conditional aligned report totals:
`DirectGeometryCandidate=32`, `PseudorealReindexCandidate=8`,
`CalibratedTripletCGCandidate=6`, `CorrectedReusedInvariantCandidate=4`,
`slots=50`, `blocks=46`. The simple unrotated geometry report becomes
`ExactGeometryCandidate=36`, `ExplicitGeometryMismatch=8`,
`TripletBasisPending=6`.

Do not trust stale pre-fix exports as post-fix evidence. Existing
`stage22_summary.json` remains *Uncertified* by design.
