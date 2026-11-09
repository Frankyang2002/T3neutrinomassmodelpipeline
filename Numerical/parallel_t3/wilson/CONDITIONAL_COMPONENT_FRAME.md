# Workstream A: symbolic RGBeta invariant diagnosis and conditional component reconstruction

## Verified cause of the probe failure

The prior `scripts/RunT3ParallelAExportRGBetaComponents.wl` passed integer slot labels to RGBeta's symbolic invariants. This is not valid component evaluation. RGBeta uses Einstein summation rules, including `del[rep, a, a] = Dim[rep]`. Thus the returned `3` and `9` in the user's A_+0/A_-2 examples are **sums**, not Kronecker-delta components. Even the prior 10 `NativeInvariantComponentsEvaluated` statuses cannot be treated as certified component data.

The actual Stage22 native declaration inventory contains the group-invariant heads `del`, `eps`, `delS2`, plus the project-defined `T3RGBetaS2TripleInvariant` (three `delS2` and three `eps` tensors) and repeated dummy triplet index `a`. These structures have definite *abstract* contraction identities, but the orientation of `eps` and basis vectors of the S2 representation remain coordinate choices. The dedicated Python interpreter evaluates the source-exported InputForm expressions in a **declared**, orthonormal highest-weight symmetric-pair basis (fundamental indices `1,2`; triplet components `11, (12+21)/sqrt(2), 22`; `eps[1,2]=+1`). It verifies the triplet orthonormality/projector identities and rejects all other AST constructions. It does not call RGBeta's native component evaluator and must never claim to have recovered RGBeta's original concrete CG phase conventions. No source files belonging to B or C are changed.

The new JSON output is **not** an input to `RunT3ParallelATensorDictionary.py` as a normalized potential: the native `Coupling` function, `UniqueArrangements` permutations, `+/- Lagrangian` signs, `PlusHc`, scalar field unitary rotations, and the anti-Hermitian A(-4) portal still require certification. The original Matchete tensor export has its own CG frame, and equal SU(2) dimensions do not imply equal component coordinates. The new file's `prefactor_certified=false` and `native_triplet_frame_certified=false` flags are hard blockers to a physical mapping.

## Run from the repository root

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
.\.venv\Scripts\python.exe -m unittest discover -s tests -p "test_t3_parallel_a*.py" -v
.\.venv\Scripts\python.exe -m scripts.RunT3ParallelAReconstructRGBetaComponents `
  --declarations output\general_t3_stage22\rgbeta_native_declarations.json `
  --registry RGE\general\generated\general_rgbeta_quartic_registry.json `
  --output output\parallel_t3\wilson\rgbeta_abstract_components_conditional_v2.json
```

If `output\general_t3_stage22\rgbeta_native_declarations.json` does not exist locally, rerun the already present `scripts\ExportGeneralT3ScalarTensorRGBetaDeclarations.wl` to generate it; do **not** substitute the old `rgbeta_native_components_probe.json`.

Expected result (subject to actual declaration content): 16 models, 50 reconstructed invariant arrays, no unsupported symbolic heads. Failure is explicit and non-destructive. No existing output is overwritten. This batch was tested with synthetic and source-excerpted declarations; the user's complete local 50-slot JSON was not available in this environment, so **16/16 actual input execution has not been verified here**.

## Next physical step

To determine a Matchete↔RGBeta dictionary, derive the transformation of *all* fields' generators and invariant tensors into a common real scalar basis, evaluate the RGBeta **fully normalized quartic action** (including `Coupling` and `UniqueArrangements`), then compare potential polynomial columns and their ranks/holdout residuals with the Matchete source. Explicitly retain both A(-4) Higgs-triplet portals: `lambdaH2Inv1` is anti-Hermitian before multiplication by `I` inside `PlusHc`, and `Inv2` is Hermitian. Neither the signed GammaCC-to-Weyl factor nor complete one-loop Wilson ADM is resolved by this reconstruction.
