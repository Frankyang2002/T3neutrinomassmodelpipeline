# T3 full-stack overnight campaign v1 (additive candidate calculation)

## Objective and strict boundary

This package **implements and connects the calculable numerical pieces** of the current T3 main-branch code without modifying Workstream A's Lagrangian or the production RGBeta/Matchete solvers. The full physically renormalized 16-model UV→EFT matching→neutrino-mass problem is **not solved** by the package. It deliberately does not fill unavailable Wilson thresholds with zero or reuse the restricted B/C one-loop triangle for a new, independently nonzero mixed bilinear.

Source inspected: `Frankyang2002/T3neutrinomassmodelpipeline`, public main commit `7ebcef87d92a` (2026-10-31 13:17:56 UTC). A subsequent local Workstream A edit may invalidate saved RGBeta exports; the preflight rejects stale source hashes. Since the complete Windows checkout and Wolfram installation are not accessible to the package builder, the source-backed full-campaign computation **must be run locally**.

## New files (complete, not patches)

| Destination relative to project root | Role |
|---|---|
| `Numerical/experimental/T3GeneralMatrixBeta.py` | New **full 3x3** extra-Yukawa state on top of the native production `StateVectorLayout`, 38+12 quartic registry, 16-model RGBeta expression evaluation, exact loop factor handling |
| `Numerical/experimental/T3BCThresholdInteractionBridge.py` | Existing mixed-doublet mass matrix / epsilon eigenbasis / SU(2) / arbitrary ordered scalar tensor and bilinear action transport checks across the B/C trajectory |
| `scripts/RunT3FullStackOvernight.py` | Campaign orchestration, 16-model source preflight, UV full-matrix ODEs, actual B/C full-mass UV/EFT1 runner, v23 source-normalized matching reference, optional historical Y/X pipeline, per-model checkpoints and convergence |
| `scripts/RunT3FullStackOvernight.ps1` | Overnight PowerShell entry point |
| `tests/test_t3_full_stack_overnight_v1.py` | Synthetic regression, negative controls and optional live-source parser/source-inventory tests |
| `docs/T3_FULL_STACK_OVERNIGHT_V1.md` | This technical contract and report interpretation |

The package adds files only. It does not edit any existing `T3RGBetaModel.wl`, Matchete Lagrangian, production RGE, matching routine, benchmark or GitHub repository.

## What each executable stage actually does

**1. All 16 coupled UV RGEs**: reads `Reports/output/workstream_a/integrated_general_yukawa_uv/native/<model>.json`, the Stage29 *candidate* quartic initial values and existing fitted-reference UV inputs. Uses the original `Numerical.core.RGBetaEvaluator.evaluate_inputform_expression`, `StateVector` packing and `BetaVector` loop-factor conversion. It evolves **all exported base beta variables plus all added quartics**; in the six eligible A/D/E cases, it evolves the independent complex `y1Conj` or `y2Conj` **3×3 matrix** in the SAME ODE. A nonzero extra-Yukawa seed (default `0.3` times an existing Yukawa magnitude) is a **test-only RGBeta-coordinate parameter**, not a physically fitted coupling. Neither the Matchete/RGBeta map nor the extra eight holomorphic scalar quadratic beta functions are certified by this first-stage export.

It also locates candidate heavy-fermion thresholds from the *undecoupled* matrix-valued UV flow, using native SVD/Takagi singular masses. These roots are **not** physical sequential-decoupling scales after the first fermion is removed.

**2. Actual B/C five-real-mass UV→EFT1**: calls `scripts.RunT3BCCoupledUVToEFT1.run_model` as currently implemented, with native UV and EFT1 full-mass beta exports; includes running `mHSq`, `mS1Sq`, `mS2Sq`, `Re(muT3GeneralS1S2)` and `Im(muT3GeneralS1S2)`, plus the original complete RGBeta dimensionless sector. This is a **separate trajectory from stage 1**: B/C do not have the extra-Yukawa matrix in this class. Fermions are removed as a grouped **tree-level renormalisable** projection at 100 TeV, with no complete Wilson matching.

**3. Eigenbasis and scalar interaction transport**: reuses `MixedDoubletMassBasis`, `quadratic_mass_tensor_transport`, `verify_tensor_transport` and `verify_su2_covariance` to check the running `X/Y` spectrum, inverse field transformation, holomorphic bilinear plus H.c., five generic source interaction field patterns and consistency through the grouped-F boundary. These are actual mathematical mass/tensor transport checks but **not** a certified cross-package native CG coupling dictionary.

**4. Restricted B/C leading fixed-order matching**: invokes the existing Workstream C v23 Matchete-based full-mass all-heavy-versus-sequential comparison, verifies current four source-audit SHA256 digests, and checks the physical-C5 flag remains false. It does **not** pretend this frozen-input reference uses the new B/C stage-2 RG background. A nonzero independent mixed quadratic (stage 2 default) is outside the archived v23 matching source. Full finite wavefunction/parameter/EOM matching is outstanding.

**5. Existing conditional EFT1/Y/X/SM numerical transport**: `--include-legacy-sequential` optionally replays the original `RunT3BCYThresholdProjection` then `RunT3BCEFT2Sequential` with new nonoverwriting output IDs. This is the **historical `continuous_v1` trajectory**, *not* the new B/C five-mass stage-2 trajectory; it still omits finite Y and X threshold terms, and therefore remains conditional. These outputs are not combined with the v23 frozen-input matching reference. No fake callback is registered in the new `GeneralT3EFTHandoffController`.

**6. Overnight validation**: strict beta inventory/source-status and source-sha preflights, source-checked resume, per-model failure isolation, scalar polynomial and SU(2) tests, and `--convergence` second UV integration at `rtol/5`. The source-backed campaign completes only after it runs on the user's full checkout; a synthetic unit test is not an empirical physics certificate.

## Installation and exact PowerShell commands

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "$env:USERPROFILE\Downloads\t3_full_stack_overnight_v1.zip" `
  -DestinationPath . -Force

python -u -X utf8 -m unittest discover -s tests -p 'test_t3_full_stack_overnight_v1.py' -v
python -u -X utf8 -m scripts.RunT3FullStackOvernight --preflight

# The long campaign: 16 coupled UV systems + B/C full mass and field maps
# + restricted matching reference + optional older Y/X conditional replay.
& .\scripts\RunT3FullStackOvernight.ps1
```

**If Workstream A modifies `RGE/running/rgbeta/T3RGBetaModel.wl` after the saved UV export**, first regenerate the 16 native RGBeta export results with the existing Workstream A stack:

```powershell
python scripts\RunT3ExtraYukawaStack.py --all
python -u -X utf8 -m scripts.RunT3FullStackOvernight --preflight
```

The *source* and beta input SHA256s must agree. Do not modify or suppress the check. A changed Stage29 registry similarly requires regeneration of Stage29 evidence before interpreting initial quartics.

For a faster single-model check:

```powershell
python -u -X utf8 -m scripts.RunT3FullStackOvernight `
  --all --models B_-1 --phase uv --output-id t3_full_stack_smoke_v1
```

For **zero optional Yukawa boundary** as a control:

```powershell
python -u -X utf8 -m scripts.RunT3FullStackOvernight `
  --all --phase uv --extra-seed-fraction 0 `
  --output-id t3_full_stack_extra_zero_v1 --convergence
```

Existing outputs cannot be overwritten silently. Use `--resume` only with identical source and numerical settings; for a changed source/seed/settings choose a new `--output-id`. A failure is recorded and later models continue unless `--fail-fast` is specified. A `progress.json` is updated after each model.

## Output locations and interpretation

- `output/studies/t3_full_stack_overnight_v1/summary.json`: global campaign status and **explicit physical limitations**
- `output/studies/t3_full_stack_overnight_v1/status.csv`: per-model stage status
- `output/studies/t3_full_stack_overnight_v1/uv/<model>_full_matrix_uv.json`: extra 3×3 Yukawa and quartic trajectories, source hashes, UV heavy-mass candidate crossings, convergence if selected
- `output/studies/t3_full_stack_overnight_v1/bc_full_mass/<model>_uv_to_eft1_coupled_five_mass_v1.json`: existing B/C live RGBeta trajectory
- `output/studies/t3_full_stack_overnight_v1/mass_basis/<model>_mixed_mass_basis.json`: running mass/tensor transport consistency checks
- `output/studies/t3_full_stack_overnight_v1/matching/bc_restricted_v23.json`: archived **frozen-input** matching and low-energy proxy reference
- `output/studies/t3_full_stack_overnight_v1/logs/`: subprocess output, including matching and optional Y/X replay
- `output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/full_stack_legacy_y_v1/` and `full_stack_legacy_x_v1/`: older **conditional** `continuous_v1` two-threshold outputs

Expected flags include `physical_C5_certified=false`, `full_16_model_sequential_matching_performed=false`, and `physical_sequential_matching_certified=false`. These are intentional, not software failures.

## What remains physically missing

1. Workstream A's final source-normalized Matchete/RGBeta UV Yukawa, mass and quartic component dictionary (and full UV Lagrangian equivalence).
2. Inclusion of *all* eight nonbaseline scalar mass bilinears and all six extra Yukawas in a universal single full-UV RG state, using the independently validated full mass beta sources. The current 16-model UV export has the extra Yukawas but explicitly lacks extra bilinears; B/C stage 2 instead includes the mixed mass bilinear.
3. Consistent EFT1/EFT2 *Wilson* evolution and complete finite field, EOM, and parameter matching across individually removed F/Y/X thresholds, for general nonzero additional operators.
4. General scalar and Yukawa physical threshold matching for the other fourteen models.
5. An independent Feynman-rule/global-sign normalization and full three-family EFT anomalous dimension for a strictly certified `C5`.

The batch **does not** claim these unsolved physics problems are fixed by running longer. Its purpose is to finish the large implementable, source-controlled numerical infrastructure and generate actionable missing-data failures, while retaining the valid restricted Matchete results.
