# General T3: independent SU(2) audit of four remaining reused B/C quartics

## Why this stage is needed

Stage 23 returned 50/50 native RGBeta prefactors, 46/46 full-rank **conditional** coupling matrices, and `physical_coupling_mapping: UNCERTIFIED`. During the next global field-basis check, no single conjugation-consistent `I` or epsilon-`J` transformation on the three complex doublets `(H, S1, S2)` aligned all six B/C quartic blocks with Matchete simultaneously. The maximum was five of six blocks. The reported **eight pseudoreal-reindexed slots** are only per-block tensor identities, not a verified global field basis.

The actual **registered** B/C RGBeta expressions (from the Stage 22 native Wolfram exports) reveal that four *other* reused coupling definitions, each used in **both** `B_-1` and `C_-1`, fail the explicit SU(2) generator test under the same fundamental/conjugate convention that verified the corrected `110101` invariant.

## Exact candidate contractions

For each `AddQuartic`, component indices `1,2,3,4` refer to the ordered fields in the declaration; `eps` is the antisymmetric SU(2) fundamental invariant with `eps[1,2]=+1`, and `del` is a Kronecker delta. The replacements below are **proposed**: they are not yet installed in `RGE/running/rgbeta/T3RGBetaModel.wl`.

| Registered coupling | RGBeta fields | Candidate `GroupInvariant` function (schematic component form) |
|---|---|---|
| `lambdaS1barS2S2bar2` | `(S1†, S2, S2†, S2†)` | `(ε₁₄δ₂₃ + ε₁₃δ₂₄)/2` |
| `lambdaS1bar2S2bar2` | `(S1†, S1†, S2†, S2†)` | `(ε₁₄ε₂₃ + ε₁₃ε₂₄)/2` |
| `lambdaS1S1bar2S2bar` | `(S1, S1†, S1†, S2†)` | `(δ₁₃ε₂₄ + δ₁₂ε₃₄)/2` |
| `lambdaHHdagS1barS1bar` | `(H†, H†, S1†, S1†)` | `(ε₁₄ε₂₃ + ε₁₃ε₂₄)/2` |

The full native Mathematica expressions are in `scripts/ExportGeneralT3BCQuarticNativeCovariance.wl`; the explicit index-by-index Python expressions are in `scripts/BCQuarticRemainderCovariance.py`.

For the two B/C models, the **Python exact calculation** establishes:

- Eight original source-extracted SU(2) tensors have nonzero Lie-generator residuals; all eight proposed tensors have zero residual for `T1`, `T2`, and `T3`.
- Each candidate tensor equals the corresponding native Matchete 1D tensor **component by component**, without an additional global `J` reindexing: `R=M`, `A=[1]`.
- Each candidate agrees with the *previous conditional* per-block `J` witness. This does not verify that `J` is a globally valid field redefinition.
- A restricted scan of all `2³ = 8` common transformations that apply `I` or `J` to both members of each doublet/conjugate pair finds **no** global alignment of all six B/C blocks with the unmodified four reused invariants. This scan is not a proof about arbitrary unitary basis transformations.
- These are geometric and covariance results only. **No physical coupling map, scalar-potential sign, or beta function is certified.**

## Evidence integrity

The Python script calls the existing `QuarticPrefactorAudit.run()` first; this checks 16 models, 50 native slots, 46 blocks, registered RGBeta prefactors, Matchete tensor ranks, and the aligned geometry SHA-256 input hashes. It then checks the exact names, fields, and source expressions of the eight remaining B/C reused slots. It **fails closed** if the native expressions are stale or have already been changed in an unreviewed way.

The independent Wolfram script builds the real RGBeta `B_-1` and `C_-1` models, loads the generated extension and reads their **registered** `$quartics` invariant functions using `RGBetaNativeAccess.wl`, then evaluates both the native and candidate contractions under the three generators. The component heads `del`/`eps` are converted before substituting numerical indices to avoid RGBeta's `del[rep,a,a] = Dim[rep]` dummy-index trace convention being misapplied to component labels.

The live exporter has **not yet been run in this environment**. After it runs, `VerifyGeneralT3BCQuarticNativeInvariants.py` checks all sixteen original/proposed component arrays (two models × four couplings × two versions) and all covariance statuses against the independent Python evaluation.

## Files in this ZIP (complete new files)

- `scripts/BCQuarticRemainderCovariance.py`
- `scripts/ExportGeneralT3BCQuarticNativeCovariance.wl`
- `scripts/VerifyGeneralT3BCQuarticNativeInvariants.py`
- `tests/python/test_general_t3_bc_covariance.py`
- `docs/GENERAL_T3_BC_REUSED_COVARIANCE.md`

**No production file, Stage 22 script, Stage 23 script, or GitHub repository is changed by this package.**

## PowerShell commands (repository root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage24_BC_Covariance.zip -DestinationPath . -Force

python -m pytest -q tests\python\test_general_t3_bc_covariance.py
python -m scripts.BCQuarticRemainderCovariance

wolframscript -file .\scripts\ExportGeneralT3BCQuarticNativeCovariance.wl `
    output/general_t3_stage22/stage24_bc_wolfram_native_covariance.json

python -m scripts.VerifyGeneralT3BCQuarticNativeInvariants
```

Expected Python result: 8/8 native failures, 8/8 covariant proposed tensors, no restricted global `I/J` field basis, physical mapping **UNCERTIFIED**. The live verifier should report `native_failed_expected: 8`, `proposed_su2_pass: 8`, `exact_python_wolfram_components_agree: True`. If a check fails, **do not** apply the candidates as production changes. Please provide the live output JSON for source review.

## Next code change, after successful live verification

Modify only the four shared `GroupInvariant` functions in the **complete local** `RGE/running/rgbeta/T3RGBetaModel.wl` replacement. Then regenerate native RGBeta declarations/invariants and update Stage 22–23 geometry and prefactor diagnostics to classify these eight slots as *direct corrected-invariant candidates*, not pseudoreal reindexings. Rerun the local tests and keep physical normalisation uncertified pending a genuinely native common real-scalar action comparison.
