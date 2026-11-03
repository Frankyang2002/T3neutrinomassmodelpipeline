# B2 G2 stage 2 — installed RGBeta symbolic generator/field-index witness

## Scope and source provenance

Source reviewed on public GitHub at `Frankyang2002/T3neutrinomassmodelpipeline` main `2ac58cc3de12c7af11e2e7e214415dadf0138f4f` (read-only). Reviewed, unchanged source blobs include:

- `RGE/running/rgbeta/T3RGBetaModel.wl`: `7af7a4ddf99a7e38210543a21b3b1cda203eaf02`
- `Lagrangian/interactions/T3Topology.wl`: `a45a099a8b766cbad0114a0559cb0bcb80895649`
- `Lagrangian/model/SU2Invariants.wl`: `dbe0cd490204dec363fcf4278f079b58c8ddcefc`

The runtime launcher reuses `source_provenance(root)`, which requires **all six** reviewed physics-source blobs to match Git HEAD, index and canonical clean-filter worktree. It passes validated HEAD and three B2 source hashes via environment to Wolfram. No Git command runs in the Wolfram child. No fixed HEAD is required: future commits are allowed only if the pinned physics-source blobs remain unchanged.

Upstream RGBeta v1.2.0 source previously inspected:

- `aethomsen/RGBeta/Package/FieldsAndCouplings.m`, Git blob `d524fa948fdbcf70ea5a61358b13a78b16a38faa`.
- `aethomsen/RGBeta/Package/GroupsAndIndices.m`, Git blob `1be60d95d7c6ef6dbde2585c5a811a5284db99d2`.

The upstream implementation writes the scalar gauge generator into conjugate-paired field blocks. For an unbarred scalar carrying `gRep`, its symbolic `Tscal` entries involve the two orientations

\[
  (\overline\phi,\phi): \quad +\tfrac12\,t^{A}_{ij}(gRep),
  \qquad
  (\phi,\overline\phi): \quad -\tfrac12\,t^{A}_{ji}(gRep).
\]

RGBeta also encodes the generator rule

\[
 t^A_{ij}(\overline{R})=-t^A_{ji}(R).
\]

These are **upstream symbolic source definitions**, not an observed numerical doublet/triplet axis. The native probe independently checks the loaded `tGen[Bar[SU2L[fund]],A,i,j]` identity and captures `DownValues` of the installed `Tscal` symbol.

## What is measured at runtime

The Python launcher starts an isolated Wolfram kernel, registers restricted A–E representatives using existing `T3RGBetaModel.wl`, and captures:

1. Installed RGBeta package path and loaded native `Tscal` definition.
2. Whether its symbolic `tGen[Bar[SU2L[fund]],A,i,j]` is exactly minus the reversed unbarred generator.
3. Scalar and fermion gauge representations `S1`, `S2`, `F`, plus their field-index-map positions. Both `S_i` and `Bar[S_i]` must have **distinct** scalar slots in the live registry.
4. Native `y1`/`y2` field orders and quartic `lambdaT3` registration in all five classes.
5. Whether `eps[SU2L[fund],1,2]` and `delS2[SU2L,1,1,2]` still remain symbolic.

The Python validator refuses wrong native status, missing/duplicate classes, changed source inventory, incorrect SU(2) dimensions, missing conjugate scalar slots, false generator identity, and all attempted physical-certification flags.

## What this does *not* resolve

Neither `Tscal` symbolic generator blocks nor `tGen[Bar[...]]` establishes a numerical realization of RGBeta's abstract invariant `eps` or `delS2`. Thus this batch cannot determine their global sign/phase or a triplet weight-axis map into Matchete. It does not convert the G1 internal factors `sqrt(2),2,6` into physical couplings. The prior exact B/C conditional intertwiner remains a convention hypothesis, not a measured RGBeta numeric tensor. B1 `GammaCC` remains separate.

**Physical certification: 0/16. Safe for production integration: NO. No EFT/RGE/one-loop matching performed.**

## Corrected native package-private index registry access

The first Windows probe (`_01`) ran its mock-based Python tests but failed in
Wolfram: `Lookup::invrl` on unevaluated `$fieldIndexMap`. Inspection of
upstream RGBeta `Package/FieldsAndCouplings.m` shows that `$scalars`,
`$fermions`, `$yukawas` and `$quartics` are `PackageExport`, whereas
`$fieldIndexMap` is **`PackageScope`**. The unqualified identifier in our
`Global` probe therefore did not reference the loaded package registry.

The corrected probe requires **exactly one** native name matching
a name beginning with `RGBeta` and ending with `$fieldIndexMap`, evaluates it only after model registration, and
requires the resulting value and its `"All"` entry to be associations.
The emitted JSON records the resolved symbol name and each family's positive,
distinct scalar / conjugate scalar registry slots. This is an observed native
package-index registry, not numeric SU(2) component indices or an action map.

The diagnostic `_01` file is immutable and is not accepted as a success.
All new output uses `_02`.

## Output

New, read-only witness: `output/general_observable/physical_normalisation/b2_rgbeta_generator_native_02/orientation_witness.json`. Its recorded source HEAD and byte hash come from the fresh Windows run. If the native source structure differs, the script writes a *diagnostic* JSON with status `BLOCKED_NATIVE_SYMBOLIC_GENERATOR_WITNESS` and exits nonzero. Inspect the diagnostic; do not infer a source bug.

## Run on Windows

```powershell
Set-Location "C:\thesistemp\T3neutrinomassmodelpipeline"
$zip = Join-Path $env:USERPROFILE "Downloads\T3_B2_G2_Native_Private_Index_Map_Fix.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
.\.venv\Scripts\python.exe -m pytest tests\test_physical_b2_rgbeta_generator_native.py -q
if ($LASTEXITCODE -ne 0) { throw "G2 Python tests failed" }
$env:T3_EXPERIMENTAL_GENERAL_FERMIONS = "0"
$env:T3_EXPERIMENTAL_GENERAL_QUADRATICS = "0"
$env:T3_EXPERIMENTAL_GENERAL_QUARTICS = "0"
$env:T3_EXPERIMENTAL_FLAVOR_MASSES = "0"
.\.venv\Scripts\python.exe scripts\RunPhysicalB2RGBetaGeneratorNative.py `
  --output "output\general_observable\physical_normalisation\b2_rgbeta_generator_native_01\orientation_witness.json"
if ($LASTEXITCODE -ne 0) { throw "G2 native symbolic generator probe failed" }
```

Do not invoke the Wolfram `.wl` directly: the Python process supplies the measured source provenance. Do not overwrite the existing output. If the `_02` path is occupied, use a fresh `_03` suffix. Never overwrite the `_01` diagnostic. Existing G1 `_05` and G2 conditional `_01` outputs remain unchanged.
