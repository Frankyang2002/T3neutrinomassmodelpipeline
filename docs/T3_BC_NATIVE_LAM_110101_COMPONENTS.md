# B(-1)/C(-1): live RGBeta `Lam` component-coefficient gate

## Why this stage exists

The previous live probe on the local machine **succeeded** after correcting the
symbol context to `RGBeta`PackageScope`Lam`. It reported `sparse arrays=1` and
both `110101` quartic couplings visible for B(-1) and C(-1). It did **not**
extract the actual coefficients of the native `Lam` tensor. This package adds
that narrow independent check without changing the original probe or any RGE.

The check reads the **installed** RGBeta v1.2.0 `Lam` tensor, takes its one
native sparse component array, and selects the ordered scalar slots
`{H, Bar[H], Bar[S1], Bar[S2]}` and the Hermitian-conjugate ordering
`{Bar[H], H, S1, S2}`. It isolates the two registered quartic couplings
`lambdaT3GeneralQ110101I1` and `lambdaHHdagS1barS2barCross`, and checks
all 16 SU(2) entries in each of the four (coupling, orientation) sectors.

RGBeta's `AddQuartic` constructs `24(1/Sqrt[2])^4=6` for four complex fields.
The stored native data records `24` distinct arrangements for the 110101
four-distinct-field case. **Hypothesis tested in this stage:** a selected
ordered entry of `Lam` equals `(6/24)=1/4` times its registered invariant,
for direct and conjugate orientations. The Wolfram script does not assume
this succeeded: it compares all entries and exports the first eight mismatches.
It exits nonzero if any component fails. The factor `1/4` is a coefficient
*inside the symmetric real-field tensor*, not the physical coefficient of
the complex-field quartic potential.

## Files

- `scripts/ProbeT3BCNativeLam110101Coefficients.wl`: new, direct native RGBeta
  `Lam` evaluation and SU(2) component comparisons.
- `scripts/VerifyT3BCNativeLam110101Coefficients.py`: new, strict verifier of
  the live result.
- `tests/test_t3_bc_native_lam_110101_coefficients.py`: new regression tests.
- `docs/T3_BC_NATIVE_LAM_110101_COMPONENTS.md`: this record.
- `docs/T3_GENERAL_SCALAR_MATCHING_RESOLUTION_LOG.md`: updated master log.

## Install and run (PowerShell, repository root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_BC_native_lam_component_extraction_v4.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/test_t3_bc_native_lam_110101_coefficients.py
wolframscript -file scripts/ProbeT3BCNativeLam110101Coefficients.wl output/audits/t3_bc_native_lam/native_lam_coeff_v4.json
if ($LASTEXITCODE -ne 0) {
  if (Test-Path output/audits/t3_bc_native_lam/native_lam_coeff_v4.json) {
    Get-Content output/audits/t3_bc_native_lam/native_lam_coeff_v4.json
  }
  throw "Native Lam component comparison failed; inspect diagnostics"
}
python -u -X utf8 -m scripts.VerifyT3BCNativeLam110101Coefficients `
  --input output/audits/t3_bc_native_lam/native_lam_coeff_v4.json `
  --output output/audits/t3_bc_native_lam/native_lam_coeff_v4_verified.json
```

A pre-existing output is never overwritten; select another `vN` filename.
This stage may expose unexpected RGBeta component convention or combinatorial
normalization. **Do not change the expected factor or relax the checks to
force a pass.** In that case, share the first mismatch and examine the native
`Lam`/`CouplingTensorMap` definitions to establish the correct projection.

## Scope and status

- Verified by this stage **only if the user-local Wolfram audit succeeds**:
  the 128 individual real-valued component comparisons in the two named
  doublet quartic sectors, using the installed RGBeta `Lam` implementation.
- Previously verified by separate stages: Matchete native extra-quartic mass
  rotations, the 2x2 `110101` geometry map, and a source-reconstructed
  Matchete-versus-RGBeta potential polynomial.
- **Not certified:** a live Matchete action comparison against the RGBeta
  `Lam` tensor, the complex-to-real field projection at the action level,
  independent vertex/Feynman-rule matching, Yukawa terms, the other
  quartic blocks, physical general scalar matching, and sequential EFT
  neutrino mass predictions.
- The physical certification flags in this package remain false. No production
  source, default benchmark, past audit output, or GitHub files are modified.

## Runtime extraction correction / diagnostic (31 October 2026)

The user's RGBeta v1.2.0 run of the original native coefficient probe reached
`T3BCSelectCoupling` but failed with `Lam component has no native Matrix terms`.
The 12 Python regression tests passed; those tests do not inspect the installed
Wolfram sparse-array contents. No `native_lam_coeff_v1.json` was written.
The earlier independent `native_lam_v2.json` diagnostic did demonstrate that
both 110101 coupling names occur somewhere in the native `Lam` tensor, but it
did **not** prove that they occupy the two *requested* sparse-array entries.

The updated probe:

- Uses a safer nested-head `Matrix` predicate and traverses expression heads
  when finding `Matrix[coupling][...]` terms.
- Records `Extract[array, positions]`, the scalar field index map, sparse-array
  dimensions, the runtime `Matrix` symbol names, and excerpts of actual
  coupling-bearing sparse-array rules. This discriminates a zero entry/index
  mismatch from a coupling-term recognition failure.
- Writes **diagnostic JSON even when a model is blocked**, then exits with a
  nonzero code. The unchanged Python verifier accepts **only** the earlier
  successfully completed 64-component certificate for **each** model. It must
  never be run as if the diagnostic failure were a pass.

The source-level factor `6/24 = 1/4` remains a **hypothesis under test** at the
native tensor level; no numerical coupling map or matching code was changed.

Run from the repository root with an unused output filename:

```powershell
python -X utf8 -m pytest -q tests/test_t3_bc_native_lam_110101_coefficients.py
wolframscript -file scripts/ProbeT3BCNativeLam110101Coefficients.wl output/audits/t3_bc_native_lam/native_lam_coeff_v2.json
if ($LASTEXITCODE -ne 0) {
  Write-Host "Native coefficient check blocked; inspect diagnostic JSON:"
  Get-Content output/audits/t3_bc_native_lam/native_lam_coeff_v2.json
  throw "Live coefficient check not certified"
}
python -u -X utf8 -m scripts.VerifyT3BCNativeLam110101Coefficients --input output/audits/t3_bc_native_lam/native_lam_coeff_v2.json --output output/audits/t3_bc_native_lam/native_lam_coeff_v2_verified.json
```

If the script remains blocked, provide the new JSON or the printed
`requested entry` and `bearing rules` diagnostics. Those values are needed to
correct the native sparse-array coordinate or symbol extraction using evidence
from the **actual installed** RGBeta runtime. Do not infer them from the
reconstructed potential.

## v3 correction: evaluated native `Lam` uses plain coupling symbols (31 October 2026)

The user ran the v2 native Wolfram diagnostic in their local RGBeta v1.2.0.
It blocked for **both** B(-1) and C(-1) on the extractor's assumption that
`Matrix[coupling][...]` persists inside the evaluated sparse-array entries.
The actual requested **direct** entry, reported identically in both models, was:

```text
(lambdaT3GeneralQ110101I1*del[SU2L[fund], t3a, t3d]*eps[SU2L[fund], t3b, t3c])/4
+ (lambdaHHdagS1barS2barCross*del[SU2L[fund], t3a, t3c]*eps[SU2L[fund], t3b, t3d])/4
```

There were 48 coupling-bearing sparse-array rules in each full tensor.
Therefore **the direct ordered symbolic coefficient already has the expected
`1/4` prefactor**. The v2 failure is an extraction-format error, not a
physical disagreement. It did **not** test the conjugate orientation or all
component indices.

The v3 probe uses Wolfram's `Coefficient[Expand[entry],coupling,1]` on the
*actual native evaluated sparse-array entry*, with explicit support for
`Bar[coupling]` or `Conjugate[coupling]` in the h.c. entry. It requires exactly
one nonzero conjugate encoding per coupling; reconstructs the complete
selected sparse entry from extracted terms, refusing extra or nonlinear
contributions; and then checks **2 couplings x 2 orientations x 16 indices =
64** components per model. It records the native coupling form in the JSON.

This does **not** evaluate Matchete live or certify a full physical mapping.

Install `T3_BC_native_lam_component_extraction_v3.zip` at the repository root,
then run with a fresh output path:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_BC_native_lam_component_extraction_v3.zip' -DestinationPath . -Force
python -X utf8 -m pytest -q tests/test_t3_bc_native_lam_110101_coefficients.py
wolframscript -file scripts/ProbeT3BCNativeLam110101Coefficients.wl output/audits/t3_bc_native_lam/native_lam_coeff_v3.json
if ($LASTEXITCODE -ne 0) {
  Get-Content output/audits/t3_bc_native_lam/native_lam_coeff_v3.json
  throw 'Native Lam coefficient component test not yet passed'
}
python -u -X utf8 -m scripts.VerifyT3BCNativeLam110101Coefficients `
  --input output/audits/t3_bc_native_lam/native_lam_coeff_v3.json `
  --output output/audits/t3_bc_native_lam/native_lam_coeff_v3_verified.json
```

The original v2 diagnostic remains available and is not overwritten.

## BC-NATIVE-LAM-05 v3 local failure and v4 correction (2026-10-31)

The user-local v3 run failed with `Part::pkspec1` / `Extract::psl1`, reporting
`slot$17285` (and `slot$18790`) as an invalid part specification. Both
models were blocked; no v3 JSON existed to verify. This was a **Wolfram
iterator syntax bug** in the probe. The final expression inside the outer
`Do` was separated from its iterator by `;` rather than `,`:

```wolfram
(* Incorrect: the iterator is part of the body; slot remains unassigned. *)
Do[body; {slot,1,2}]

(* Correct: the iterator is a separate argument. *)
Do[body, {slot,1,2}]
```

The complete replacement v4 script changes the separator to a comma and adds
a regression test for exactly that failure. It retains the live `Coefficient`
method, dual orientation check, 64 components/model and fail-closed
physical-certification flags. The script has not been run under installed
RGBeta by the package author; it must be executed locally with a **new** output
filename. Do not claim normalization certification until the live run passes.
