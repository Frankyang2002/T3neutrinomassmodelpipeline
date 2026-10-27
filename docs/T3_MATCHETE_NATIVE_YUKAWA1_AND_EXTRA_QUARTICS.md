# B_-1 / C_-1 native Matchete y1 and additional-quartic recoupling

Experimental mass-eigenstate Matchete source replacement, 2026-10-27.

## Verified live input motivating this change

The user's preceding run passed the diagonal free Lagrangian and the first eleven disjoint baseline interaction groups. `BaselinePart12 y1|CG:T3Y1CG` returned `$Failed` for both B_-1 and C_-1. The complete mass-eigenstate `CheckLagrangian` therefore remains unverified.

## Code changes

- **Rebuild `y1` as a native `PlusHc` Yukawa vertex.** The heavy scalar is S1 (not S2) before transport, so the scalar doublet index needs no epsilon: S1 -> U11 X + U12 Y. The vertex retains the original `NCM[Bar[l], CConj[F]]` structure and carries a single direct `T3MENativeY1CG` invariant. One independent SU(2) invariant for dF=1 and dF=3; exact original/new component equality is checked.
- **Rebuild all seven additional scalar quartic couplings with native four-leg CG tensors.** Original species are ordered H, Hdag, S1, S1dag, S2, S2dag. For each generated operator, use the original, boson-symmetrized physical tensor from `T3ExtendedBosonicTensorBasis`, contract its S2/S2dag tensor legs with canonical epsilon, and verify its nonzero rank lies in `PhysicalSU2InvariantBasis` for the new mass-field conjugation orientations. Register one direct CG per source coupling; expand X/Y mixtures only on the external field factors.
- **Preserve source inventory.** Losslessly separate both y1 and all seven extra-quartic source terms, including Hermitian conjugates. Reassemble the eigenstate Lagrangian from the existing direct-CG baseline terms, native y1/y2 vertices and native quartics, with no old scalar labels and with the same source coefficients.
- **Keep failure isolation.** Isolated baseline y1 and each of the seven extra quartics now use the exact same new expression as the full Lagrangian. Do not interpret an isolated check as a substitute for full heavy-mass-basis validation.
- **Fail-closed output schema.** Two new certificates, `native_yukawa1_tensor` and `native_additional_quartic_tensors`, are checked in Python, including all seven exact expected quartic keys, coupling labels, invariant ranks and epsilon-leg counts. Physical Matchete-to-RGBeta coupling normalization remains uncertified.

## Numerical basis and physical interpretation

For Phi=(S1,epsilon conjugate(S2))^T = U (X,Y)^T, both X and Y are heavy, complex SU(2) doublets with Y=-1/2. Source-leg transformations are:

```
S1       = U11 X + U12 Y
S1dag    = conj(U11) Xdag + conj(U12) Ydag
S2       = -eps [conj(U21) Xdag + conj(U22) Ydag]
S2dag    = -eps [U21 X + U22 Y]
```

The SU(2) rank/epsilon checks are mathematical group-tensor tests only; these do not certify the physical phase/normalization between RGBeta and Matchete.

## Verification performed

- Local Python regression: **111 passed** (includes the unchanged mixed-scalar numerical tests).
- Python source compile: passed.
- WL source lexical bracket/comment/string balance: passed.
- Temporary installed Matchete `InvariantTensors` check: B_-1 and C_-1 y1 each have one invariant, and for all six independent *structures* (seven couplings) the projected four-leg tensors lie in the mass-basis invariant span with rank 2, exactly. The Matchete kernel has no working front end for `CheckLagrangian`; a complete live test has **not** passed here.

## Install from project root (PowerShell)

```powershell
Expand-Archive -LiteralPath "C:\path\to\T3_Matchete_Native_Y1_And_Extra_Quartics.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
```

The next run must report both ordinary `CheckLagrangian=True` and detailed `HeavyMassBasis=True` to pass. No threshold Wilson coefficients are computed and neither numerical EFT1 data nor numerical mass rotation needs regeneration.
