# Batch 2: independent one-loop scalar hard stage (ordinary models)

## Purpose

This package adds an **independent**, source-backed native Matchete scalar-only
one-loop matching attempt for `D_+2` and `A_-4`. The first fermion threshold is
**never rerun**: the stage imports the already exported full-general EFT1
`eft1_after_F_tree_inputform.txt` with both heavy scalar species retained.

The stage uses the general UV interaction registry, restores the matching CG
seeds, removes the integrated fermion from the field registry, promotes only
`MS1` and `MS2` from light to heavy EFT power counting, checks the heavy scalar
quadratic sector against `FreeLag`, then calls `CheckLagrangian` followed by
`Match[..., EFTOrder -> 5, LoopOrder -> 1]`. If successful it processes the
independent hard contribution using the same Matchete order as D(-2):

1. `GreensSimplify`
2. `EOMSimplify`
3. `EvaluateLoopFunctions`
4. `ReplaceEffectiveCouplings`

The output is a *hard-sector C5 candidate*. It does NOT combine the inherited
one-loop EFT1 operators, the direct running log, or the conditional LLSS scalar
bubble from Batch 1. It has NOT been tested with Mathematica in this environment
and must not be described as completed until its local reports show success.

Both supported models are **ordinary scalar-mass** cases with no extra quadratic,
quartic, or conjugate Yukawa interactions in the full-general inventory. The
mixed models B(-1)/C(-1) require physical mixed-scalar eigenfields; the six
self-bilinear models require a real mass basis. No ordinary-mass shortcut is
used for those cases.

## Files and intended project paths

- `Lagrangian/RunT3GeneralOrdinaryScalarHard.wl`
- `scripts/RunT3GeneralOrdinaryScalarHard.py`
- `tests/python/test_t3_general_ordinary_scalar_hard.py`
- `README_T3_SEQUENTIAL_BATCH2.md`

The ZIP contains **complete files**. None replace the successful D(-2) pipeline
or the previous Batch 1 numerical runner. Do not install from the ZIP inside a
pre-existing source folder: extract into the **repository root**.

## Local PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_Sequential_Batch2_OrdinaryHard.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_t3_general_ordinary_scalar_hard.py
python -X utf8 -m scripts.RunT3GeneralOrdinaryScalarHard --preflight-only
```

If both saved source boundaries are ready, run the first model only:

```powershell
python -X utf8 -m scripts.RunT3GeneralOrdinaryScalarHard `
  --models D_+2 --run-id ordinary_hard_01
```

Or both ordinary models, without repeating the F matching:

```powershell
python -X utf8 -m scripts.RunT3GeneralOrdinaryScalarHard `
  --models D_+2,A_-4 --run-id ordinary_hard_02
```

To reuse any successful model while retrying a blocked one:

```powershell
python -X utf8 -m scripts.RunT3GeneralOrdinaryScalarHard `
  --models D_+2,A_-4 --run-id ordinary_hard_02 --resume
```

The runner writes only to
`output/general_sequential_pilot/ordinary_scalar_hard/<run-id>/`.
For each model, inspect `hard_scalar_input_audit.json`,
`hard_scalar_report.json`, `hard_c5_candidate.txt`,
`wolfram_stdout.log`, and `wolfram_stderr.log`.

**Success marker:** `hard_scalar_report.json` must have
`Status = HardScalarCandidateExported`. A successful hard candidate is still
not the full `C5`: inherited + running pieces and NF=1 pole cancellation
must be assembled and checked separately. The absolute SU(2) spinor
normalization remains provisional.

## Inputs and non-destructive behavior

The ordinary models use saved full-general F-only sources from:

- `output/general_observable/fermion_first_all16_v1/D_p2/data/`
- `output/general_observable/fermion_first_all16_v1/A_m4/data/`

The Python preflight requires source metadata identifying all fermions as
integrated and both scalars as still active, a 3-flavor F boundary and both
PL/PR LLSS chiral sectors. It hashes all six key source files and this WL
backend, refuses unrecognized operator additions, and prevents overwriting a
successful hard match after source inputs change.

Do not replace these newer full-general sources with the older one-scale
`output/general_matchete_matching/` results or with an earlier restricted
`output/full16_symbolic/` threshold without re-validating the interaction
inventory and matching conventions.

## Validation status

- Python preflight/unit tests: provided and runnable without Wolfram.
- Wolfram one-loop Match: **not run in this packaging environment**.
- D(+2), A(-4) independent scalar hard matching: **not yet verified**.
- B(-1) mixed scalar hard matching: **not included**.
- Full sequential C5: **not yet assembled** for these models.


## Repair for ordinary_hard_01 failures (2026-11-07)

The local logs showed two distinct blockers:

- `A_-4`: `wolframscript` treated `-4` as an option and returned usage exit 2.
  The Python launcher now passes `m4`, which the Wolfram script parses as -4.
- `D_+2`: `CheckLagrangian` emitted `FindDummyIndices::trippleindex` for
  `y1[p,r] y2[q,r] / MF[r]` in the F-first tree (three occurrences of the
  same NFlavor index). The new Wolfram stage leaves the input tree unchanged
  but changes the matching *working copy* into a one-heavy-mass-eigenvalue
  symbolic template by replacing indexed `MF[r]` by symbolic `MF[]`, with
  an exact count audit. It then calls `RelabelIndices` to avoid collisions
  between dummy names used by Flavor and NFlavor.

This yields a **per-heavy-eigenvalue hard scalar matching template**, not a
full three-mass flavor C5. The numerical lift must evaluate the template at
all distinct heavy masses with the corresponding `y1[:,r]` and `y2[:,r]`
and sum those contributions. The already successful conditional LLSS run
must not be overwritten or added to a one-scale C5 as an extra term.

Extract this repair ZIP over the project root, then run:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m pytest -q tests/python/test_t3_general_ordinary_scalar_hard.py
python -X utf8 -m scripts.RunT3GeneralOrdinaryScalarHard --preflight-only
python -X utf8 -m scripts.RunT3GeneralOrdinaryScalarHard `
  --models D_+2,A_-4 --run-id ordinary_hard_01 --resume
```

Previously successful hard outputs are only reused with matching file hashes;
failed models are retried. New Wolfram `CheckLagrangian` and `Match` success
cannot be established without running the scripts on Windows.
