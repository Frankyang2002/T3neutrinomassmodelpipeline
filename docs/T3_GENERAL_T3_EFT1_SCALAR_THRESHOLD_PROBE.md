# Ordinary T3 `B_-1`/`C_-1`: scalar-threshold Matchete probe

## Purpose and status

The **already live-passed** first threshold, `F` heavy with eigen-doublets `X,Y` light, produces a *tree-level* EFT1 interaction and (in a separate run) an unrenormalized one-loop hard `LLHH` sector. The hard-pole audit found 68 finite, 68 matching-log and 68 pole terms for each model. Matching those sectors coefficientwise to each other is an internal audit, **not** a hard-pole / EFT1 RGE test.

This batch tests the *second* field-elimination boundary, independently and without treating the first hard contribution as renormalized:

    EFT1(tree: X,Y light, F absent)
       -> fresh Matchete process with F absent and X,Y heavy
       -> Match(EFTOrder=5,LoopOrder=0) and Match(EFTOrder=5,LoopOrder=1)
       -> raw one-loop scalar-threshold LLHH sector

This is **not** the full fermion-first sequential EFT calculation. In particular, it does **not** propagate the one-loop hard contribution, compute EFT1 Wilson tensor running, remove UV divergences through independently checked EFT1 counterterms, split X/Y matching thresholds, or report physical `C5`.

## Actual implementation

- The one-loop first-threshold JSON and matching files are validated using the existing first-threshold certificate code; the prior heavy-mass-basis and tree matching certificates are also required.
- The **tree-only** `*_eft1_tree.inputform.txt` is taken from the successful fermion-first one-loop export (not the `tree_plus_loop` input), preventing a second one-loop counting of inherited hard terms.
- A new `wolframscript` process loads Matchete, the original model and the same native SU(2) CG construction used in the earlier live-passed scripts. It removes `F` from the field registry, preserves its mass as a real coupling, and re-registers `X,Y` as `Mass -> {Heavy,MSX/MSY}`.
- The saved tree EFT1 is parsed in the fresh process and requires ordinary and detailed `CheckLagrangian`, including `HeavyMassBasis=True`.
- Tree and one-loop matching are run on **the same tree EFT1**. Their difference is exported alongside the matched SM tree+loop EFT and its raw `LLHH` sector.
- The raw selection requires two `H` and two `l` field entries and the Matchete `GammaCC` bilinear, and excludes `X`, `Y` and `F`. A zero raw sector is reported as zero, not silently promoted to physical matching.
- Progress stages and logs survive an aborted/failed Wolfram run. Existing Lagrangian sources, saved EFT1 outputs, RGEs and numerical benchmarks are not edited.

## Project paths

New source files, installed at repository root:

- `Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl`
- `scripts/RunGeneralT3EFT1ScalarThresholdProbe.py`
- `tests/python/test_general_t3_eft1_scalar_threshold_probe.py`
- `docs/T3_GENERAL_T3_EFT1_SCALAR_THRESHOLD_PROBE.md`

Output directory:

`output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe/`

Per model it contains `*.json`, `*.progress.json`, `*.log`, `*_sm_tree.inputform.txt`, `*_sm_tree_plus_loop.inputform.txt`, `*_scalar_one_loop_correction.inputform.txt`, and `*_raw_scalar_llhh.inputform.txt`.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_General_EFT1_Scalar_Threshold_Probe.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_eft1_scalar_threshold_probe.py
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model B_-1
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model C_-1
```

This needs the successful previous `RunGeneralT3FermionFirstEFT1Tree` and `RunGeneralT3FermionFirstEFT1OneLoop` results in their original output locations, the same mass-eigenstate input parameters, and your working local Matchete + `wolframscript`. A Matchete failure is expected to be diagnosed using `last_stage` in the `.progress.json` and the final `.log`; do not claim success from offline unit tests.

## Next stage if the probe passes

Build an **independent mass-eigenfield EFT1 tensor map** so that the exact Matchete `LLSS` and transformed scalar quartics populate the intermediate Wilson/renormalisable RGEs in the **same SU(2) basis**; then test whether the direct `LLSS -> LLHH` mixing and the hard-pole residue satisfy the specified MSbar consistency condition. Incorporate its RG evolution and the required counterterms before assembling the physical symmetric `C5` at the scalar threshold and running to 100 GeV.
