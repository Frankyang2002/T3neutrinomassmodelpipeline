# T3 EFT1 scalar-threshold: verified renormalizable reconstruction and quartic isolation

## Status

This is a **fail-closed experimental repair** of `SavedEFT1TreeHeavyScalarCheckFailed`, not a demonstrated physical threshold-matching solution. It has not run against local Wolfram/Matchete in this environment. Do not infer successful two-threshold C5 until the user's Wolfram runs succeed and the subsequent renormalization/RGE checks are complete.

The local B/C logs show `Hermiticity=False` in the saved EOMEFT even when X/Y are still light, while canonical normalization, heavy-mass basis, charge/gauge checks and field inventory pass. The earlier sector audit implicates combinations containing scalar self-quartics but does not isolate an individual coupling. Neither changing the scalar classification nor `GreensSimplify` / `HcSimplify` cures it.

## Verified reconstruction attempt

`Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl` rebuilds the *renormalizable* SM + X,Y interactions from the same native SU(2) mass-eigenstate components as the first-stage UV Lagrangian. It keeps the **exact** saved 20 (`B_-1`) or 24 (`C_-1`) tree-matched Majorana LLSS terms, and removes the F Yukawas/propagating field. It constructs the candidate:

`L_rebuilt = L_SM + FreeLag[X,Y] + L_rotated_scalar_interactions + L_LLSS_saved`.

Crucially, **this candidate is not used unless** `(i)` `Expand[L_EOMEFT - L_rebuilt]` vanishes exactly, or strict symbolic simplification after index relabelling establishes equality to zero; `(ii)` *ordinary* `CheckLagrangian` returns `True`; `(iii)` the full detailed Matchete `Hermiticity`, `CanonicallyNormalized`, `HeavyMassBasis`, `GaugeInvariant` and `AllObjectsDefined` checks all return `True`. A further detailed check is performed after acceptance. No `PlusHc` is added, no terms are removed or made zero for matching, and no projected `(L+HC[L])/2` is accepted.

If the candidate is algebraically different from the saved EFT—even because of a legitimate yet unproven EOM identity—the runner **refuses** it. This may make the repair conservative and unsuccessful; a no-go result is diagnostic, not evidence the EFT is physically non-Hermitian.

The exact report is `<stem>_reconstruction_audit.json`; if accepted the success certificate and `<stem>_hermiticity_diagnostics.json` record `VerifiedRenormalizableReconstruction`. The Python runner verifies the proof flags before accepting the success certificate.

## Quartic coupling isolation (diagnostic only)

On failure, the runner partitions all terms exactly as before and additionally groups the 122 scalar self-quartic terms by Matchete `Coupling[lambda..., ...]` label. For each label it checks `background + coupling` and `background + all self-quartics except coupling` and writes detailed flags under `scalar_self_quartic_coupling_isolation` in `<stem>_hermiticity_sector_audit.json`. These temporary diagnostic deletions **never** feed `RunT3Matching` and cannot change a physical coefficient.

## Installation and execution (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_General_EFT1_Verified_Reconstruction_Fix.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_eft1_scalar_threshold_probe.py
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model B_-1
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model C_-1
```

Rerunning fermion-first tree/one-loop, prior hard-pole audit and earlier numerical results is unnecessary.

## If Matchete still fails

Inspect the two generated files:

- `output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe/B_m1_final_eft1_scalar_threshold_probe_reconstruction_audit.json`
- `output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe/B_m1_final_eft1_scalar_threshold_probe_hermiticity_sector_audit.json`

and corresponding `C_m1` files. In particular look for `algebraic_equality_proven`, `independent_matchete_check_passed` and `scalar_self_quartic_coupling_isolation`. These distinguish a true operator mismatch from a Matchete representation issue and identify the coupling to correct in the next source-level change.

## Scope

Even if the second scalar-threshold matching runs, its exported LLHH sector remains unrenormalized and cannot yet be combined with the fermion threshold without EFT1 running, matching-scale/pole checks and proper flavour normalization. The SM-only one-scale neutrino benchmark remains the currently verified numerical endpoint.
