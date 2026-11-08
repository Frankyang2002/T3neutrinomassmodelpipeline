# B/C(-1) X-only Matchete post-export report recovery v4

## Evidence actually inspected

The user-supplied `xy_fresh_native_evidence_v3.zip` contains 148 records and the exact old Wolfram source with SHA256:

`ccf16e4c5aac78a8fc80c5f50fac766bec65b62279c40ae99db97cbe8e698ee0`

B and C each show:

- All six `CheckLagrangian` conditions true for the heavy-free-rebased X-only matching input.
- Both `RunT3Matching[matchingInput,5,0]` and `[matchingInput,5,1]` traverse all five native stages and log `Matching pipeline completed successfully.`
- Checkpoints include `x-followup-one-loop-matching-returned`, `x-followup-operator-export-started`, `x-followup-matching-exports-returned` and `x-fresh-X-followup-returned`.
- Saved nonempty `x_tree_sm`, `x_full_sm`, `x_loop_difference`, `x_raw_llhh` InputForm files (36 independent additive LLHH terms each). These are **raw pole-containing expressions**, not renormalized finite coefficients.
- The only observed post-export status-update failure is `AssociateTo::argrx`; the old source called `AssociateTo[report,rule1,rule2,rule3,rule4,rule5]` instead of using `Join[report,<|rules|>]`.
- The old report remains `XFollowupOneLoopMatchingInProgress` even though the Matchete stages and operator exports completed; this is a reporting defect, not an X Match failure.

## Changes

`Lagrangian/RunT3BCXFromSavedY.experimental.wl`: full replacement of the exact uploaded 119,547-byte v3 source. The three incorrect multi-rule `AssociateTo` calls (tree-failure branch, loop-failure branch, success branch) are corrected to `report=Join[report,<|...|>]`. The final status export now fails closed. All mass-basis construction, quartics, SU(2) CGs, EFT2 field assignments and calls to Matchete are unchanged.

`scripts/RecoverT3BCXRawMatching.py`: NEW read-only evidence recovery. It validates the original v3 source hash, old report, heavy X check, running log, ordered stage checkpoints and all seven saved operator files. It counts the 36 raw LLHH additive terms and creates new recovery reports and an evidence ZIP **without modifying the old v3 results and without rerunning Wolfram**. The recovered status is independent of the broken v3 status; no false physical certification flags are raised.

`tests/test_t3_bc_x_raw_recovery_v4.py`: new tests including missing-file, wrong source hash, omitted one-loop checkpoint, bad matching-input fallback, and two-model atomicity negative controls.

## Installation on the user's local repository

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\t3_bc_x_raw_match_recovery_v4.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p "test_t3_bc_x_raw_recovery_v4.py" -v
python -u -X utf8 -m scripts.RecoverT3BCXRawMatching --all
```

Expected:

- `B_-1: X tree=PASS; X one-loop=PASS; raw LLHH terms=36; physical=NOT CERTIFIED`
- `C_-1: X tree=PASS; X one-loop=PASS; raw LLHH terms=36; physical=NOT CERTIFIED`

Outputs at:

`output/audits/t3_bc_x_raw_matching_recovery_v4/`

- `B_m1_x_raw_recovery.json`, `C_m1_x_raw_recovery.json`
- `summary.json`, `summary.csv`
- `recovered_x_native_raw_evidence_v4.zip`

**Do not rerun Matchete** merely to fix the v3 status; the raw outputs are already present. The Wolfram replacement fixes subsequent native matching experiments.

## Physical limitations

1. The F-first saved EOMEFT is not algebraically certified equivalent to the independently reconstructed Y-only EFT1.
2. The X matching input is the tree-level Y EOMEFT. The finite Y one-loop correction is not inherited, so full sequential matching has not been performed.
3. Source parameter substitutions at running scales `muY` and `muX` have not been made.
4. The raw X LLHH expressions contain dimensionally regulated poles `ϵ`. UV subtraction, matching between EFTs, pole/log scale consistency and finite coefficient extraction remain necessary.
5. The conditional numerical Wilson running retains existing one-generation and scalar mapping limitations. No physical neutrino spectrum is certified.

**Next substantive physics batch:** use the recovered native Y and X raw expressions to do a pole-and-finite operator decomposition, compare UV/IR pole structures with the appropriate EFT counterterms, and only then map a renormalized C5 boundary onto the `continuous_v1` numerical running trajectory. Do not conflate frozen source parameters with running scalar boundaries.
