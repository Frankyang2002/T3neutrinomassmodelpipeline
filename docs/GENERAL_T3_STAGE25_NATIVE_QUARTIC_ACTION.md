# Stage 25 — Native RGBeta complex-quartic assembly probe (read-only)

## Scope and result status

Stage 24 repaired the eight remaining reused B/C invariant slots and its live
native post-fix audit reported eight exact matches. Stage 23 gave 50/50 native
quartic prefactors and 46 *conditional* signed maps. **None of these is a
certification of the physical Matchete–RGBeta potential mapping.**

Stage 25A checks an independent part of the missing normalisation step: the
**actual RGBeta runtime quartic tensor assembly**, not merely the stored
`Coupling` and `UniqueArrangements` strings. It also checks RGBeta's complex
scalar species index map and scalar contraction metric. This is a necessary
cross-check, but the exporter intentionally does **not** claim to have obtained
a Cartesian real-scalar component tensor or established a potential sign.

## Source conventions checked

Checked against the upstream `aethomsen/RGBeta` `Package/FieldsAndCouplings.m`
source (current repository at preparation):

- `UpdateFieldIndexMap` (around lines 130–155) represents each complex scalar
  `phi` and `Bar[phi]` as separate index entries, with scalar contraction
  metric 2 for the conjugate pair.
- `AddQuartic` (around lines 751–787) constructs
  `normalization=24 (1/sqrt(2))^4=6` for four complex slots and calls
  `CouplingPermutations`.
- `Package/GroupsAndIndices.m` (lines 275–277) contains
  `Matrix[m_][] = m`.  Thus the stored function
  `6*Matrix[lambda][##1] &` evaluates at zero flavour indices to the
  symbolic coupling **`6*lambda`**, not the unevaluated
  `6*Matrix[lambda][]`.  Stage 25 checks this exact expression in Wolfram
  and independently checks its serialisation in Python.  This change fixes
  only the runtime-format mismatch; it does not modify physical normalisation.
- `AvgPermutations` / `CouplingTensorMap` (around lines 925–943)
  populate the tensor using each unique field-index permutation with
  weight `1/N`.
- `Lam` (around lines 1035–1046) builds its quartic tensor by using
  `CouplingTensorMap` and explicitly including the conjugate term when
  `SelfConjugate=False`.

The live exporter uses the *runtime* `CouplingTensorMap`, not a reimplementation
substituted in its place. For each of the 50 Stage 22 native declarations, it:

1. Rebuilds the relevant RGBeta model and loads the actual generated extra
   quartic declarations (16 models).
2. Confirms that the registered `Fields`, `Invariant`, `Coupling`, and
   `UniqueArrangements` match the refreshed Stage 22 declarations.
3. Evaluates `CouplingTensorMap` and independently reconstructs its ordered
   rules from the native permutation list; compares each LHS and full symbolic
   RHS (including field-index permutation and SU(2) tensor factors).
4. Checks that all four scalar slots have a conjugate partner in the native
   field-index map and that each pair has metric component 2.
5. Checks the exact native complex-quartic coupling factor.

The Python verifier additionally recalculates Stage 23 (including the Stage 22
source-file SHA256 consistency checks), compares it to the saved Stage 23 report,
and checks all 50 live permutation counts across the 46 conditional blocks.

## What this does NOT show

`$fieldIndexMap` uses complex field species and conjugate species as independent
bookkeeping entries; it does **not** directly label the Cartesian real and
imaginary components. In particular, the factor `6` in `Lam` alone does not
prove the sign or physical coefficient in the Matchete Lagrangian.

The following are still uncertified:

1. A single globally consistent Matchete–RGBeta SU(2) component intertwiner,
   particularly the independent triplet basis and relative phases.
2. The physical action sign connecting a Matchete `+lambda * operator` in
   `LUV` to RGBeta's real/quartic action convention. Stage 23's
   `lambda_R=-A^{-1} lambda_M` remains an **assumption-dependent candidate**.
3. Direct comparison of the **entire native quartic action**, with the
   Hermitian-conjugate contribution and identical-field combinatorics, after
   converting complex scalar species to an independently verified common
   Cartesian basis. Stage 25A tests RGBeta's assembly engine but not the
   `Lam` contraction against a separately exported physical Matchete action.
4. The eight previously blocked scalar quadratic entries.

No production RGE, vertex definition, numerical benchmark, or certification
flag is modified by this package.

## Installation (PowerShell from repository root)

Place `T3_General_Scalar_Stage25_NativeAssembly.zip` in
`C:\thesistemp\T3neutrinomassmodelpipeline` and run:

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage25_NativeAssembly.zip `
    -DestinationPath . -Force

python -m pytest -q tests\python\test_general_t3_stage25_quartic_action.py

wolframscript -file .\scripts\ExportGeneralT3Stage25NativeQuarticAssembly.wl `
    output/general_t3_stage22/stage25_live_native_quartic_assembly.json

python -m scripts.Stage25QuarticActionAudit
```

Prerequisites: regenerated Stage 22 native declarations and aligned geometry,
Stage 23 prefactor report after Stage 24 corrections; `RGBeta` available in
`wolframscript`. If these have not been run, regenerate before invoking Stage 25:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage22.ps1 `
    -ExportRGBetaDeclarations -NativeDiagnosticsOnly
python -m scripts.Stage22RGBetaAlignedGeometry
python -m scripts.Stage23QuarticPrefactorAudit
```

Successful Stage 25 live run writes:

- `output/general_t3_stage22/stage25_live_native_quartic_assembly.json`
- `output/general_t3_stage22/stage25_native_quartic_assembly_diagnostics.json`

Expected verified coverage: **16 models, 50 slots, 46 blocks**. The total
number of native permutation rules is measured at runtime; it is intentionally
not predicted here.

If the Wolfram run fails, retain the original error output and share it;
do not interpret the lack of a JSON report as a successful certification.
The Wolfram exporter was not executable in the preparation environment.

## Next numerical decision

If the live assembly passes, the remaining blocking substep is to compare a
minimal set of independently oriented, Hermitian-conjugated quartic *actions*
with the same component basis and fix the relative sign. Apply the mapping to
the RGE only if that explicit comparison is successful. For rough benchmarks,
limit this check to the couplings actually used at the target point rather
than expanding to unnecessary phenomenology.
