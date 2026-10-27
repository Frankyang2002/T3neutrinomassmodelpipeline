# Native SU(2) mass-basis tensors for S2-dependent baseline scalar quartics

## Scope

Experimental `B_-1` and `C_-1` native Matchete mass-eigenbasis validation.
This package replaces only the following full project files:

- `Lagrangian/RunGeneralT3MassEigenMatchete.wl`
- `scripts/RunGeneralT3MassEigenMatchete.py`
- `tests/python/test_general_t3_mass_eigen_matchete.py`

This document is new at `docs/T3_MATCHETE_NATIVE_S2_SCALAR_QUARTICS.md`.
No GitHub files are changed. The existing numerical mass-basis JSON and EFT1 output are unchanged.

## Live failure that motivated the change

Both models passed the native reconstructed `lambdaT3` and `lambdaH2Inv1/Inv2` tests but failed `DefineCG::notInvariant` for `BaselinePart8 lambda12Inv1|CG:T3CrossScalarCGInv1`. The old four-leg CG multiplied by independent epsilon bridge CGs is not always accepted by the final Matchete invariant check.

## Implementation

The original source Lagrangian is strictly partitioned into the already-native topology and H2 portals; all remaining `lambda12Inv1`, `lambda12Inv2`, and `lambdaS2Inv1` terms are then separated with an exact algebraic reconstruction check. No source coupling is deleted without native reconstruction.

For `lambda12Inv1/Inv2`, the old tensor with legs `(S1bar, S1, S2bar, S2)` is contracted with epsilon in its last two slots to obtain the direct mass-basis invariant with legs `(Xbar, X, X, Xbar)` (where `X` denotes either eigenfield `NewScalarX` or `NewScalarY`). The two projected tensors and the new invariant basis must both have rank 2, and their union must also have rank 2.

For `lambdaS2Inv1`, all four old scalar legs are epsilon-transported from `(S2bar,S2,S2bar,S2)` to `(X,Xbar,X,Xbar)`. Both source and target are projected onto their correct symmetric-boson subspaces; the single projected tensor must lie in the single native invariant span.

The exact symbolic coefficients follow `S1 = U11 X + U12 Y`, `S2bar = - U21 epsilon X - U22 epsilon Y`, and their conjugates. For all three quartics the two/four minus signs cancel. No Matchete/RGBeta normalization map is inferred.

The actual new CG tensors are registered using `DefineCG` and the new basis's `CGReps`. The fail-closed result JSON now includes `native_cross_scalar_epsilon_tensors` and `native_S2_self_epsilon_tensor`. Python validates their counts and exact-rank flags alongside the existing topology and H2 certificates. Isolation diagnostics use the same reconstructed expressions as the final Lagrangian.

## Verified and not verified

- The Python regression suite passed (99 tests in this package's combined local suite).
- Wolfram-Language delimiters and Python syntax compiled.
- A temporary native Matchete kernel returned projected/native span ranks `(2,2)` and `(1,1)` and exact Bose symmetry for the self tensor.
- **Not verified:** full `CheckLagrangian` and `HeavyMassBasis` in the user's Matchete instance; no numerical threshold or Wilson matching has been performed.

## Run

From the repository root (PowerShell):

```powershell
Expand-Archive -LiteralPath "C:\path\to\T3_Matchete_Native_S2_Scalar_Quartics_Fix.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
```

For a failure, the runner prints the first failing isolated operator. The corresponding full log is in `output/general_t3_mass_eigen_matchete`.
