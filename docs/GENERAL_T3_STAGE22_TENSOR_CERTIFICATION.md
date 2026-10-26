# Stage 22: RGBeta native association-key repair (Fix 2)

## Current status and verified evidence

**Physical Matchete–RGBeta quartic normalisation remains UNCERTIFIED.**

Four JSON files from the user's Windows run were inspected (no source files were modified):

- `matchete_tensors.json` contains exactly 50 explicit Matchete SU(2) tensors in 16 models. The independent Python verifier passes 46 full-rank field-content blocks, including the four two-invariant blocks. There are 26 triplet-containing slots, and six tensors contain explicit imaginary components; those must not be mapped by assuming identical triplet or conjugate component conventions.
- `rgbeta_native_declarations.json` has 50 coupling *names*, but the actual `Fields`, `Invariant`, `Coupling`, `UniqueArrangements` values are all serialized as `Missing["KeyAbsent", ...]`.
- `rgbeta_native_invariants.json` repeats those unresolved values; `Missing["KeyAbsent", Invariant][...]` is **not** a native invariant evaluation.
- `stage22_summary.json` consequently gives misleading positive diagnostic labels: the original validators only checked that metadata strings were nonempty. These cannot be used as coupling certification evidence.

The root cause is a **Wolfram symbol-context mismatch**: RGBeta associates native metadata under package-context symbol keys. The previous exporter looked them up with unqualified `Fields`, `Invariant`, `Coupling`, etc., which may resolve to different `Global`` symbols outside the package. This diagnosis is supported by the missing-key payload and RGBeta source construction, but the repaired Wolfram exports still require execution on the user's Windows installation.

## Fix 2 implementation

Full replacement files, all paths relative to repository root:

1. `scripts/Stage22RGBetaNativeAccess.wl` — **new**; resolves native `$quartics` Association keys by their symbol name, checks exactly one match, insists on native four-field lists, and rejects `Missing` values before serialization.
2. `scripts/ExportGeneralT3Stage22RGBetaDeclarations.wl` — **replacement**; uses the common accessor, fails closed and clears stale output before a new export.
3. `scripts/ExportGeneralT3Stage22RGBetaInvariants.wl` — **replacement**; evaluates the actual native stored invariant function, rejects missing/unresolved expressions, and clears stale output.
4. `scripts/RunGeneralT3Stage22.py` — **replacement**; rejects serialized Wolfram `Missing[...]` in both channels; validates explicit native four-field lists, compares corresponding native fields/invariants/coupling functions/arrangements between files, and marks all 50 slots uncertified on malformed evidence.
5. `tests/python/test_general_t3_stage22.py` — **replacement**; regression coverage for the exact failure seen in uploaded JSON and cross-file metadata mismatch.
6. `docs/GENERAL_T3_STAGE22_TENSOR_CERTIFICATION.md` — **replacement**; this document.

No change to native Matchete exporter, scalar operator generation, RGBeta model, production numerical trajectories, mass matrices, or threshold code.

## Installation and tests (PowerShell, project root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage22_Fix2.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_stage22.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage22.ps1 `
  -ExportRGBetaDeclarations -NativeDiagnosticsOnly
```

This deliberately **reuses the already validated** `output/general_t3_stage22/matchete_tensors.json`. It reruns the two repaired RGBeta diagnostics and checks every field-content block. `-ExportRGBetaDeclarations` in the Stage 22 PowerShell runner also exports native invariant expressions.

Expected verified diagnostics after a successful native run:

- `output/general_t3_stage22/matchete_tensors.json`
- `output/general_t3_stage22/rgbeta_native_declarations.json`
- `output/general_t3_stage22/rgbeta_native_invariants.json`
- `output/general_t3_stage22/stage22_summary.json`

The JSON summary must contain:

- `native_matchete.status = MatcheteNativeComponentsAndRanksVerified`
- `native_rgbeta_declarations.status = NativeRGBetaDeclarationsVerifiedNoCommonBasisComponents`
- `native_rgbeta_invariant_expressions.status = RGBetaNativeInvariantExpressionsVerifiedNotComponents`
- `native_rgbeta_coherence.status = TwoNativeRGBetaDiagnosticChannelsAgree`
- `certification_status = Uncertified`, `normalisation_certified = false`

If a required native key cannot be resolved, the exporters exit nonzero, print the available key names, and do not claim success. They require actual Mathematica/RGBeta and were not run in the container (WolframScript unavailable). Local Python tests: **16 passed**.

The strict command without `-NativeDiagnosticsOnly` still fails until `output/general_t3_stage22/rgbeta_native_tensors.json` contains source-verified explicit component tensors with field alignment, potential-prefactor/conjugation provenance, and an invertible exact tensor map. **Do not** fill that file with candidate invariants or assert proof flags merely to pass tests.

## Next physical task

After obtaining the corrected native RGBeta invariant functions and field order, derive explicit RGBeta components for `del`, `eps`, `delS2`, and triplet structures, verify the intertwiners against Matchete's complex triplet basis, and establish RGBeta's potential prefactors (including the internal `AddQuartic` symmetry normalisation). Then test all 46 blocks componentwise, including B/C `110101` and E's two-invariant sectors. Only a genuine component/prefactor proof may activate the 50 transformations in the numerical RG pipeline. Eight quadratic mixing slots remain separate work.
