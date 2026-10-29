# General T3: four shared B/C RGBeta quartics — local correction

**Scope:** repair only four verified, native SU(2)-noncovariant `AddQuartic` invariant definitions in `RGE/running/rgbeta/T3RGBetaModel.wl`. This file is the full replacement of current GitHub `main` blob SHA `1711bd0c8e639347e11494013bd5b9863df0fa11`, plus precisely four expression changes. No GitHub edit.

## Inputs and evidence

The user's Stage24 local run reported **13 Python tests passed**, with **8/8 native B/C covariance failures** and **8/8 proposed SU(2)-covariant exact Matchete matches**. An independent live RGBeta export confirmed 2 models, 8 expected native failures, 8 proposed passes, and exact Wolfram/Python component agreement. This tests the standard `U`/`Conjugate[U]` field representations, not signed physical coupling normalisation.

## Exact replacements

Each is written for its existing registered field order; each `SelfConjugate -> False` stays unchanged. Expressions below use `D(a,b)=del[SU2L@fund,#a,#b]`, `E(a,b)=eps[SU2L@fund,#a,#b]`.

- `lambdaS1barS2S2bar2`, `{Bar[S1], S2, Bar[S2], Bar[S2]}`: `(E(1,4) D(2,3) + E(1,3) D(2,4))/2`.
- `lambdaS1bar2S2bar2`, `{Bar[S1], Bar[S1], Bar[S2], Bar[S2]}`: `(E(1,4) E(2,3) + E(1,3) E(2,4))/2`.
- `lambdaS1S1bar2S2bar`, `{S1, Bar[S1], Bar[S1], Bar[S2]}`: `(D(1,3) E(2,4) + D(1,2) E(3,4))/2`.
- `lambdaHHdagS1barS1bar`, `{Bar[H], Bar[H], Bar[S1], Bar[S1]}`: `(E(1,4) E(2,3) + E(1,3) E(2,4))/2`.

No other `AddQuartic`, coupling name, or potential prefactor is altered. The earlier `lambdaHHdagS1barS2barCross` correction is preserved.

## Diagnostic compatibility

`scripts/RGBetaAlignedGeometry.py` now recognises both the historical `J`-reindexed tensor candidates and the *fully* repaired Stage24 geometry profile. For the repaired profile, each of these eight slots has an **exact direct tensor equality** `R=M`, without a per-block `J`.

`scripts/QuarticPrefactorAudit.py` accepts the old and repaired **whole-model** classification profiles. It still emits a provisional signed matrix and retains `signed_physical_coupling_map_certified = false`. The shared potential-sign assumption, live real-scalar tensor normalisation and global field convention still require independent verification.

The old `BCQuarticRemainderCovariance.py`, `ExportGeneralT3BCQuarticNativeCovariance.wl`, and `VerifyGeneralT3BCQuarticNativeInvariants.py` are historical **pre-fix** diagnostics; do not rerun them against repaired native exports. Use the new post-fix exporter and verifier below.

## Windows PowerShell (run from project root)

Make a backup and install all complete replacement/new files:

```powershell
Copy-Item .\RGE\running\rgbeta\T3RGBetaModel.wl .\RGE\running\rgbeta\T3RGBetaModel.wl.pre-stage24.bak
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage24_BC_LocalFix.zip -DestinationPath . -Force
```

Check the Python regression suite:

```powershell
python -m pytest -q tests\python\test_general_t3_bc_postfix.py tests\python\test_general_t3_prefactors.py tests\python\test_general_t3_tensor_aligned_geometry.py
```

Recheck the **live registered** quartics and compare all 8 SU(2) tensors with Matchete:

```powershell
wolframscript -file .\scripts\ExportGeneralT3BCQuarticCovarianceCorrection.wl `
    output/general_t3_stage22/stage24_bc_postfix_wolfram_native_covariance.json
python -m scripts.VerifyGeneralT3BCQuarticCovarianceCorrection
```

Refresh the complete Stage22/23 native diagnostic evidence, in this order:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3ScalarTensorCertification.ps1 -ExportRGBetaDeclarations -NativeDiagnosticsOnly
python -m scripts.RGBetaAlignedGeometry
python -m scripts.QuarticPrefactorAudit
```

Expected Stage22 **conditional** geometry classes after the full correction:

```text
DirectGeometryCandidate:                 32
CorrectedBCReusedInvariantCandidate:      8
CalibratedTripletCGCandidate:             6
CorrectedReusedInvariantCandidate:        4
slots:                                   50
blocks:                                  46
```

The Stage23 signed coupling-map count remains 46. Neither the live SU(2) verification nor these exact component matches certify the physical Matchete–RGBeta coupling map; all certification flags remain false.

## Synthetic checks performed before packaging

- Source block fingerprint: only the four previously `delS2 * delS2` contractions replaced; previous 110101 repair preserved.
- Exact SU(2) test: proposed tensors pass all three generators.
- Synthetic Stage22 post-fix registry fixture: 32 + 8 + 6 + 4 = 50 slots, 46 blocks.
- Synthetic Stage23 prefactor fixture: 50/50 native factor checks, 46 conditional signed matrices.
- Python post-fix verifier regression: accepts exact in-memory component tensors and rejects corrupted SU(2) residuals, incorrect components, wrong Matchete tensors and old-format exports.

**No live post-fix RGBeta run was performed in the packaging environment.**
