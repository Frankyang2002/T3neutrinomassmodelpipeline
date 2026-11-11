# Worker C — Matchete scalar-threshold source audit (2026-11-12)

## Install scope

Additive files, **do not overwrite** any previous Worker C module or solver:

| File | Project-relative location | Role |
|---|---|---|
| Source audit | `Numerical/joint_scalar_matching/worker_c_matchete_audit.py` | Read-only evidence and mandatory physics gates |
| Runner | `scripts/RunT3WorkerCScalarAudit.py` | CLI for existing Matchete profile files |
| Regression tests | `tests/test_t3_worker_c_matchete_audit.py` | 5 synthetic offline tests |
| This report | `WORKER_C_MATCHING_AUDIT_README.md` | Verified findings, boundaries, usage |

Nothing changes in `Numerical/joint_scalar_matching/worker_c_joint.py`, the original `core.py`, the fermion threshold, Worker B code, or production orchestration. Neither GitHub nor historical outputs are modified.

## Verified from public GitHub (read-only)

Reference repository: `Frankyang2002/T3neutrinomassmodelpipeline`, `main` Git tree SHA `802659cb8e71371bebf04eff8936b4c5ae41246c` inspected 2026-11-12. **The local Windows working files were not mounted or revalidated**.

- `Lagrangian/RunGeneralDMinus2SequentialScalarThreshold.wl`: consumes the saved physical first-threshold tree and performs `Match[..., EFTOrder -> 5, LoopOrder -> 1]`, then Green/EOM/evaluated loop function/effective-coupling reduction. The script says it exports a **tentative** hard coefficient, not assembled physical SMEFT C5.
- `output/general_sequential_pilot/D_m2_stage1/data/general_d_m2_scalar_threshold/profile_tree/profile_full_scalar_one_loop_result.json`: `OneLoopRawMatchReturned=true`, `HardC5Extracted=true`, `WeinbergTermCount=68`, `FullFlavorHardC5Verified=false`, `FixedOrderSMEFTAssembled=false`, `FullGeneralScalarThresholdCertified=false`.
- `.../profile_fixed_order_components_result.json`: 34 holomorphic + 34 conjugate candidate terms, `HardC5HasY2=true` and `HardC5HasY2Conj=true`; inherited pure-SM pieces and direct running remain separate.
- `.../profile_full_scalar_one_loop_audit.json`: full saved EFT1 has 79 terms, including 8 LLSS terms: 4 `Y1Y2`, 4 `Y1Y2Conj`. Scalars registered jointly and the saved mass-power-counting audit passes; physical mass relations remain un-imposed.
- `.../profile_full_scalar_one_loop_hard_c5_candidate.txt`: explicit native `hbar` terms, poles with `ϵ`, and both `y2` and `y2Conj` couplings.
- `.../profile_c5_pole_cancellation_result.json`: after inherited assembly, simple pole cancels exactly **for `NF=1`**; `RenormalizationCertified=false`. Do not confuse this with a proof of complete renormalized matching.
- `.../profile_c5_flavour_candidate_result.json` and `.../profile_c5_normalization_reference_result.json`: flavour symmetrization and Warsaw normalization **references** exist, but absolute physical normalization remains uncertified.
- `Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl`: experimental full scalar threshold for `B_-1/C_-1` works from a tree EOMEFT and explicitly says no physical finite C5 certification. No corresponding successful saved result was located in the inspected public Git tree.
- `Numerical/parallel_t3/threshold/scalar_cw_batch3.py` computes a Coleman–Weinberg **Higgs-potential** correction, explicitly **not** C5 matching. Do not count it as an independent Weinberg contribution.

## Accounting decision

The existing `worker_c_joint.match_joint_scalars` already includes the complete finite `C_LLSS × lambda_HHSS` scalar-pair contribution within its stated truncation, plus inherited `C5_EFT1(muS)` added once. The saved Matchete scalar one-loop coefficient has **not** been decomposed into a finite `matched hard minus existing bubble` remainder with normalized flavour and the required `Y1,Y2` restriction. Hence it is incorrect to insert its raw value as `additional_scalar_one_loop`.

The audit's `safe_to_add_as_independent_c5_correction` is **always false** in this source-only module: even if future source reports show all prerequisite flags, a separate verified **numerical** remainder would still be needed. This is intentional, not a placeholder physical contribution. The module reports exactly why the source evidence is blocked and records SHA-256 fingerprints of inspected local files.

## Windows PowerShell

From the repository root with the project's Python environment activated:

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive -LiteralPath "$HOME\Downloads\T3_worker_C_Matchete_audit_20261112.zip" -DestinationPath . -Force

python -m unittest discover -s tests -p 'test_t3_worker_c_matchete_audit.py' -v

# Read-only scan of the saved D(-2) first-threshold-profile matching evidence:
python -m scripts.RunT3WorkerCScalarAudit --model D_-2

# Find and inspect the newly created report folder:
Get-ChildItem .\output\audits -Directory -Filter 'worker_c_scalar_D_-2_*' |
    Sort-Object LastWriteTime -Descending | Select-Object -First 1
```

Default profile directory, read-only:

`output/general_sequential_pilot/D_m2_stage1/data/general_d_m2_scalar_threshold/profile_tree/`

For another model, provide a compatible pre-existing profile directory:

```powershell
python -m scripts.RunT3WorkerCScalarAudit --model B_-1 `
    --profile-dir 'output\path\to\a\B_m1\one_loop_profile'
```

If that profile is missing or its schema differs, the audit records **BLOCKED** and missing/unsupported evidence rather than constructing artificial C5 terms. To reject blocked output with a failing exit code, append `--require-ready`.

A new directory is generated automatically under `output/audits/worker_c_scalar_<model>_<timestamp>_<id>/` on each invocation. The script refuses to overwrite an existing explicitly named output directory.

## Next implementation gate

For the restricted `Y1,Y2` EFT1:

1. Generate a fresh Matchete scalar-threshold result with the conjugate-scalar Yukawa *interaction itself* absent, or establish an exact physical equivalent of restricting the saved expression.
2. Isolate the **renormalized** SMEFT Weinberg coefficient including the proper EFT1 inherited pieces and soft subtraction at `muS`; determine whether any term survives subtraction of the already implemented `LLSS × HHSS` bubble.
3. Only if a nonzero independent remainder is established, export a signed, physically normalized, finite complex symmetric `(3,3)` numerical coefficient and attach its source proof. Add it to Worker C **once**, without reintroducing the EFT1 running log.

No new independent non-bubble C5 term has yet been numerically established by this audit.
