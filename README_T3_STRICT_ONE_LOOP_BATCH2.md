# Batch 2 — all-16 source-channel one-loop bookkeeping and strict EFT1→SMEFT calculation

**Project root:** `C:\thesistemp\T3neutrinomassmodelpipeline`.

**Pinned read-only GitHub main:** `4bd1230256701ac3801af4e83e8c73e5d448ee1e` (2026-11-09). **No changes to GitHub.** This package only adds new files. It does not overwrite the earlier 16/16 `source_audited_all16_v1` or `fhard_triangle_all16_batch1_v1` results.

## Physics and exact numerical scope

One **common fermion threshold** `MF=100000 GeV`: integrate out all three heavy fermions, retain both BSM scalars in EFT1. One **common scalar threshold** `MS=1000 GeV`: integrate out both scalar multiplets, then evolve the SMEFT Weinberg operator with the existing one-loop SM+Weinberg solver down to 173 GeV.

The newly implemented **strict one-loop full-Wilson endpoint** is

\[
C^{[0+1]}(M_S) = C_{LLSS}^{(0)}(M_F) + C_{5,F\text{-triangle}}^{(1)}(M_F)
 + {\ln(M_S/M_F)\over16\pi^2}\,\beta_{\mathrm{EFT1}}[C^{(0)}(M_F);g(M_F)].
\]

The numerical EFT1 anomalous dimension is the **existing** `fonseca_eq22_six_term` candidate, with all seven separately exported channels (one crossed channel is explicitly zero). It is evaluated on the zero-order **full Weyl45/real-scalar Wilson tensor**, not the one-loop hard boundary; therefore the known hard one-loop matching is not spuriously renormalized a second time at first order.

The **existing** Workstream-C `strict_one_loop` mode then takes:

- the strict full-Wilson tensor to project the Weinberg coefficient at `MS`;
- the zero-order `LLSS(MF)` tensor for the simultaneous **two-scalar eigenstate-pair bubble**;
- the exact, source-invariant **saved RGBeta potential channels** for the available benchmark quartic subset at `MS`.

This mode uses the common `MS` mass and quartic background from the saved one-loop running data, whose running represents higher-order-improved input in the scalar bubble; Workstream C itself labels that qualification. Below `MS` its existing one-loop SM+Weinberg ODE is RGE-improved; that is a separate continuation from the strict `MF→MS` expansion. Thus do not call the *entire* low-energy calculation globally strict at fixed `O(hbar)`.

**New actual calculation:** independently saves, for each model:

1. Each named source **HHSS quartic channel** contribution to the hard fermion triangle at MF, using the full saved three-heavy-fermion source-resolved tensors.
2. Each named source quartic channel contribution to the common-MS scalar bubble, by driving the existing C matcher with `quartic_real_channels_ms`.
3. Each named candidate EFT1 anomalous-dimension channel contribution to the Weinberg operator.
4. The strict first-order matched SMEFT `C5`, full SM+Weinberg RGE, masses, splittings, PMNS diagnostics, and comparisons with resummed Batch 1.

**What remains uncomputed:** all heavy-fermion one-loop hard topologies other than the quartic triangle; finite scalar-matching diagrams beyond the scalar pair bubble; scalar field normalisation/EOM matching; a certified general EFT1 ADM and signed native-to-Weyl conventions. Not physical mass predictions. Only the existing selected/reduced quartic potential (whose omissions are recorded) is decomposed; absent nonzero native quartics are **not** set to inferred physical zeros.

## Install (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\t3_strict_one_loop_batch2_all16.zip" -DestinationPath . -Force
.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_strict_one_loop_batch2.py
```

## Execute full 16-model batch

```powershell
# First full run creates a NEW output; no previous outputs overwritten.
.\.venv\Scripts\python.exe scripts\RunT3StrictOneLoopBatch2.py

# If interrupted or some models fail:
.\.venv\Scripts\python.exe scripts\RunT3StrictOneLoopBatch2.py --resume

# Print saved summary:
.\.venv\Scripts\python.exe scripts\RunT3StrictOneLoopBatch2.py --show
```

To check the first model rapidly before a full 16-model run:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3StrictOneLoopBatch2.py --model A_-4 --preflight-only
.\.venv\Scripts\python.exe scripts\RunT3StrictOneLoopBatch2.py --model A_-4 --resume
.\.venv\Scripts\python.exe scripts\RunT3StrictOneLoopBatch2.py --resume
```

`--overwrite` recomputes all selected models in **only Batch 2's own** output folder. The original B, source-audited v2, and Batch 1 outputs remain untouched.

## Output locations

- `output/parallel_t3/strict_one_loop_all16_batch2_v1/summary.json`: all models, numerical result status, missing contributions, scales, source hashes.
- `output/parallel_t3/strict_one_loop_all16_batch2_v1/strict_vs_resummed_all16.csv`: masses and mass shifts relative to Batch 1.
- `output/parallel_t3/strict_one_loop_all16_batch2_v1/inputs/<model>/strict_wilson_channels.npz`: explicit `C0`, hard boundary, strict full Wilson, and seven EFT1 `C5` running terms.
- `output/parallel_t3/strict_one_loop_all16_batch2_v1/inputs/<model>/MF_hard_triangle_by_quartic.npz`: named hard MF triangle contributions.
- `output/parallel_t3/strict_one_loop_all16_batch2_v1/inputs/<model>/quartic_channel_provenance.json`: named quartic inventories at both scales and source omissions.
- `output/parallel_t3/strict_one_loop_all16_batch2_v1/results/<model>/C5_smeft_at_MS.npz`: existing Workstream-C strict mode output, with hard + Weinberg mixing + scalar bubble.
- `output/parallel_t3/strict_one_loop_all16_batch2_v1/results/<model>/C5_per_quartic_channel.npz`: existing Workstream-C per-quartic scalar bubble results.
- `output/parallel_t3/strict_one_loop_all16_batch2_v1/results/<model>/C5_one_loop_channel_ledger.json`: explicit per-channel full complex 3x3 matrices for F-hard, EFT1 beta, scalar-MS bubble and matched `C5`.
- `output/parallel_t3/strict_one_loop_all16_batch2_v1/results/<model>/observables.json`: mass eigenvalues, mass-squared differences, PMNS where numerically defined.

No files are written to GitHub. The output hash of each source NPZ is verified before reuse. Do not delete Batch 1 outputs when evaluating Batch 2.

## Full replacement-file destinations in this ZIP

1. `Numerical/parallel_t3/wilson/strict_one_loop_batch2.py` — pure strict-first-order full Wilson endpoint and seven-channel running.
2. `Numerical/parallel_t3/wilson/quartic_channels_batch2.py` — exact selected potential channel decomposition at MF/MS and per-channel hard matching.
3. `scripts/RunT3StrictOneLoopBatch2.py` — reproducible all-16 run through existing strict Workstream C and SMEFT.
4. `scripts/RunT3StrictOneLoopBatch2.ps1` — simple PowerShell test/run wrapper.
5. `tests/test_t3_strict_one_loop_batch2.py` — mathematical and source-inventory tests.
6. `README_T3_STRICT_ONE_LOOP_BATCH2.md` — this file.

**Validation scope:** Python syntax, direct pure arithmetic test, ZIP CRC only in the generating container. The 16-model Windows numerical test has **not** been run in this container. Results from `RunT3StrictOneLoopBatch2.py` are required before asserting 16/16 completion.
