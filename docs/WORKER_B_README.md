# Worker B scalar-background handover (source-injected interface)

## Intended paths

- `Numerical/running/T3GeneralScalarBackground.py` (new; does not replace the existing modules)
- `tests/test_t3_general_scalar_background.py` (new)
- `docs/WORKER_B_README.md` (new)

## Verified upstream conventions

Read-only public reference, commit `802659cb8e71371bebf04eff8936b4c5ae41246c`:
- `Numerical/running/T3ScalarQuarticLinearResponse.py`: cached affine 12^4 real quartic, currently A(-4)/D(-2).
- `Numerical/running/T3CoupledQuarticContext.py`: source HHSS overrides candidate HHSS; physical normalization not certified.
- `RGE/running/intermediate/ScalarOnlyTensorAdapters.py`: 4 Higgs real coordinates then S1 and S2 interleaved real/imag scalar coordinates.
- `Numerical/running/GeneralT3SequentialLeadingLog.py`: 16 models with dimensions A(1,3), B/C(2,2), D(3,1), E(3,3).
- `Lagrangian/GeneralScalarMassEigenbasis.wl`: holomorphic bilinear convention V=mSq S* S + kappa S^2 + h.c. (NO 1/2), real mass matrix [[mSq+2Re(kappa),-2Im(kappa)],[-2Im(kappa),mSq-2Re(kappa)]].
- `Numerical/parallel_t3/complete/two_model_scalar_coordinate_batch16.py`: reference-to-native quartic sector maps may be incomplete.

## Interface

`construct_affine(reference, names, complex_names, direct, rank)` constructs a rank-4 quartic or rank-2 mass response from a **local** complete evaluator; checks an independent multi-coupling point.

`build_scalar_background(model, scale_gev, state, quartic, quadratic, hhss_source=..., expected_quartic_names=..., expected_quadratic_names=..., basis_rotation=...)` returns:
- `q_fixed` : (12,12,12,12) full quartic derivative tensor
- `m2_fixed`: (12,12) scalar mass squared matrix
- `rotation_original_from_mass`: (12,12), original real coordinates = rotation @ mass coordinates
- `eigenvalues_gev2`: (8,) BSM squared masses
- `q_mass`: (12,12,12,12) rotated quartic tensor
- `assumptions`: explicit assumptions/provenance

`scalar_mass_matrix` and `scalar_rotation` properties expose (8,8) arrays.

**Important:** no full 16-model local source invariant inventory, local full direct evaluator, special bilinear injection, B/C dedicated transformation, or local running state was present in this worker environment. This delivery is therefore an executable *strict integration interface*, **not** completed source-derived physics for all 16 models. It intentionally rejects unprovided inventory and B/C transformations. Do not label synthetic unit tests as benchmark validation. The caller must build both complete response bases from current local source files; do not substitute identity rotations as physical B/C input.

The EFT sequence remains UV -> all fermions out -> EFT1 with BOTH scalars -> both scalars jointly out -> SMEFT. Nothing here changes the solver or matches scalars.

## Install and test (PowerShell)

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive -LiteralPath "$HOME\Downloads\T3_Worker_B_Scalar_Background.zip" -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p 'test_t3_general_scalar_background.py' -v
```

**Next necessary source-dependent step:** construct the complete per-model quartic and quadratic `direct` evaluators using the current local general scalar invariant inventory, including D(-2) extra quadratic direction, and wire the dedicated B(-1)/C(-1) basis maps. Source-inject them into this API; then test at `mu_F=1e5` and `mu_S=1e3` GeV against the existing cached A(-4)/D(-2) results before Worker D integrates it.
