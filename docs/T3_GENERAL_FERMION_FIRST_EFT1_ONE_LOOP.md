# General T3 — first fermion-threshold one-loop matching probe

Scope: `B_-1` and `C_-1`, native scalar eigenstates `NewScalarX`, `NewScalarY`.

## Prerequisite live results

The previously supplied **tree** first-stage runner passed locally:

| Model | Matchete tree-level EFT1 `GammaCC`/`LLSS` term count |
| --- | ---: |
| `B_-1` | 20 |
| `C_-1` | 24 |

The existing `RunGeneralT3MassEigenMatchete` detailed HeavyMassBasis check and `RunGeneralT3FermionFirstEFT1Tree` report remain prerequisites. They are never overwritten.

## What this new probe actually does

The new `RunGeneralT3FermionFirstEFT1OneLoop.wl` is a **standalone full replacement-style script** based on the live-validated eigenfield and tree-first matching source. It:

1. Rebuilds the UV Matchete model with `NewFermion` **Heavy**, and both `NewScalarX`, `NewScalarY` **Light**. Both scalar eigenstates remain in EFT1. It reuses native SU(2) contraction reconstructions for both Yukawas and all 17 scalar quartics.
2. Validates the fermion-first UV Lagrangian with ordinary and detailed `CheckLagrangian`.
3. Calls the repository's `RunT3Matching[lag, 5, 0]` and validates the tree EFT1 `LLSS` output.
4. Calls `RunT3Matching[lag, 5, 1]` **on the identical UV tree-level Lagrangian**; this uses the existing Matchete stages: `Match`, `GreensSimplify`, `EOMSimplify`, `EvaluateLoopFunctions`, `ReplaceEffectiveCouplings`.
5. Computes exactly `Expand[oneLoopMatchedEFT - matchedEFT]`. The difference is the first-threshold one-loop EFT correction, avoiding inadvertent double-counting of the tree-level EFT. Checks that the correction is nonzero, explicitly includes `hbar`, and contains no `NewFermion` field.
6. Exports the tree EFT, tree `LLSS` sector, full matched EFT through one loop, subtracted loop correction, and the **raw, unrenormalized** candidate `GammaCC`/`HHLL` sector. Candidate selection requires two explicit Higgs legs, `GammaCC`, and no remaining `X`, `Y` or `F` fields. This is a conservative textual diagnostic, **not** a final operator projection.

The candidate count may legitimately be zero. That is **not** silently promoted into a physical neutrino-mass coefficient: zero is saved and reported, and interpretation requires examining the full EFT expression.

## What this does not do

The output **does not** subtract UV poles, enforce an MSbar scheme/scale, certify the hard-pole versus EFT1 anomalous-dimension consistency, establish the complete flavour/CG normalization, run EFT1 between the heavy fermion and scalar thresholds, perform scalar-threshold matching, or calculate a physical `C5`.

The correct future hierarchical computation requires, schematically,

\[
C_5(M_S) = C_{5,\mathrm{hard}}^{\overline{\mathrm{MS}}}(M_F)
 + \Delta C_{5,\mathrm{direct}}^{\mathrm{EFT1}}(M_F\to M_S)
 + \text{scalar-threshold terms}
\]

with fixed-order loop counting and the repository's Majorana symmetry convention. The expression above is only schematic; no unjustified equality between the all-heavy one-scale result and the sequential EFT result is assumed.

## Install

Extract `T3_General_Fermion_First_EFT1_OneLoop_Probe.zip` **into the repository root**, retaining the `Lagrangian/`, `scripts/`, `tests/python/`, `docs/` directories. There are no modifications to the previously successful validators.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_General_Fermion_First_EFT1_OneLoop_Probe.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_fermion_first_eft1_one_loop.py
python -m scripts.RunGeneralT3FermionFirstEFT1OneLoop --model B_-1
python -m scripts.RunGeneralT3FermionFirstEFT1OneLoop --model C_-1
```

The run uses local tree-stage and HeavyMassBasis certificates from previous successful runs. Use `--all` only after `B_-1` has completed. The new one-loop run may take longer than tree matching; `--timeout-sec` controls the external process timeout (default 7200 seconds), without modifying an existing session.

## Outputs

`output/general_t3_mass_eigen_matchete/eft1_fermion_first_oneloop/`:

- `<stem>.json`: status, stages, source and matching flags, operator counts.
- `<stem>.progress.json` and `<stem>.log`: last checkpoint and Wolfram stdout/stderr even on failure.
- `<stem>_eft1_tree.inputform.txt`: matched tree-level EFT1.
- `<stem>_llss_tree.inputform.txt`: tree-level `GammaCC` scalar sector.
- `<stem>_eft1_tree_plus_loop.inputform.txt`: matched EFT1 through one loop.
- `<stem>_eft1_one_loop_correction.inputform.txt`: exact difference between the two matched EFTs.
- `<stem>_raw_hard_llhh_sector.inputform.txt`: raw `GammaCC`/two-Higgs candidate sector; may be `0`.
- `summary_final.json`: batch summary.

## Verification status

**Offline tested:** Python runner validation, dependency checking, output schema, nonzero/missing loop factor, zero/nonzero candidate sector, status flags, mocked subprocess; combined test collection 203/203 passed.

**Not yet verified:** Actual local one-loop Matchete matching for the F-first `B_-1` and `C_-1` models. The remote Wolfram kernel cannot run the front-end-dependent `CheckLagrangian`/`Match` stages of this project's Matchete workflow.

The upstream GitHub repository is read-only for this work. Install all changes locally; do not commit or push from this batch.
