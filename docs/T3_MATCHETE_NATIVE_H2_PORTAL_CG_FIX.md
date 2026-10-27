# Native mass-eigenstate Matchete: both S2 Higgs-portal CG tensors

## Live failure (27 October 2026)

Both `B_-1` and `C_-1` pass the isolated diagonal free Lagrangian, two S1 Higgs portals, the reconstructed `lambdaT3` topology and the S1 self-interaction. The first remaining failed partition is:

```
DefineCG::notInvariant: The tensor is not invariant under group transformations.
BaselinePart6 lambdaH2Inv1|CG:T3HiggsPortal2CGInv1 -> $Aborted
```

The previous direct reconstruction of `lambdaT3` **passed the live isolated Matchete test** and is retained here.

## Change

The S2 portal includes two S2 legs. Naively transforming each field leaves the old portal CG multiplied by two epsilon bridges in the resulting Matchete expression. Both independent portal structures are now reconstructed as native four-leg CGs via the exact tensor relation

```
(Hbar, H, Bar[S2], S2) -> (Hbar, H, X_or_Y, Bar[X_or_Y])
T_new[i,j,k,l] = sum_{a,b} T_old[i,j,a,b] eps[a,k] eps[b,l]
```

The old/new basis and orientation are obtained from the existing `PhysicalSU2InvariantBasis` implementation, not from guessed Clebsch coefficients. Both original and new invariant spaces have dimension 2. The runner rejects any missing rank or span. It registers `T3MENativeH2CGInv1/Inv2`, retains the real source couplings `lambdaH2Inv1/Inv2`, and uses the existing four symbolic U entries and their conjugates. It checks exact source-term reconstruction before excluding the old portal terms and adding the projected native terms. The isolated-operator diagnostic uses the same direct native portal expressions.

The success JSON now includes `native_H2_portal_epsilon_tensors` and the Python validator rejects missing or inconsistent invariant-count, rank, and exact-projection reports.

## Verification / nonclaims

- 92/92 offline Python tests passed.
- Actual Matchete `InvariantTensors` checked in a remote Wolfram kernel: old and new invariant space dimensions both 2, projected rank 2, projected tensors lie in the native two-dimensional span.
- Full native `CheckLagrangian` is **NOT VERIFIED**: the remote kernel cannot run the front-end-dependent Matchete validation. This is a proposed fix for the user's observed live failure, not a certified live pass.
- The transformation can reveal further CG problems with other S2-dependent operators. Fail-closed diagnostics remain active.
- Numerical mass eigenvalues and U, EFT1, threshold matching and Wilson coefficients are unchanged. No physical coupling normalization is claimed.

## Full replacement project paths

- `Lagrangian/RunGeneralT3MassEigenMatchete.wl`
- `scripts/RunGeneralT3MassEigenMatchete.py`
- `tests/python/test_general_t3_mass_eigen_matchete.py`
- `docs/T3_MATCHETE_NATIVE_H2_PORTAL_CG_FIX.md` (new)

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_Matchete_Native_H2_Portal_CG_Fix.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
```
