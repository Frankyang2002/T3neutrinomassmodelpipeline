# B/C(-1) native `110101` tensor-basis map after scalar mixing

## Purpose and source provenance

This read-only audit combines **three independent saved source artifacts**:

1. `output/general_t3_stage22/matchete_tensors.json`: native Matchete `110101` tensor pair;
2. `output/general_t3_stage22/stage22_110101_wolfram_native_covariance.json` and
   `rgbeta_native_invariants.json`: native RGBeta registered source expressions,
   component tensors and live SU(2) covariance tests;
3. `output/general_t3_mass_basis/B_m1_final_mass_basis.json` and `C_m1_final_mass_basis.json`:
   existing scalar rotation matrices.

It also **requires** the previously installed
`scripts/AuditT3BCNativeQuarticMassBasis.py` for native source validation and
exact component-tensor rotation. All inputs are read only. The scripts were
inspected on the repository's main branch as of 2026-10-31. No Matchete or
RGBeta package is run during this audit; live native *saved* evidence is reused.

## Result being tested

For ordered fields `(H, Hdag, S1dag, S2dag)` and epsilon_12=+1, the two
RGBeta *invariant tensors* are

\[
R_1{}_{ijkl}=\delta_{il}\epsilon_{jk}, \qquad
R_2{}_{ijkl}=\delta_{ik}\epsilon_{jl}.
\]

`M1` and `M2` are the corresponding independent saved Matchete native
quartic tensors. The script verifies, with **exact SymPy algebra over the 16
components**, that

\[
\begin{pmatrix} R_1 & R_2\end{pmatrix}
=\begin{pmatrix} M_1 & M_2\end{pmatrix}
\begin{pmatrix}\frac12&1\\-\frac{\sqrt3}{2}&0\end{pmatrix}.
\]

The determinant is \(\sqrt3/2\), so the two independent contractions are
retained. This tensor identity was partially covered by the historical Stage22
source diagnostic; **the new check additionally tests its covariance under the
*actual evolved B/C mixed-scalar mass-basis rotation***. Each RGBeta tensor is
independently transported with the epsilon-dependent scalar field map, and
compared with the same mapped linear combination of the two Matchete tensors,
component by component and as a complex scalar polynomial.

The coupling relation `c_M = A c_R` applies **only to invariant-polynomial
coefficients in one common ordered-field and same-sign convention**. It does
**not** establish that raw RGBeta couplings are `c_R`: the native RGBeta
potential/vertex normalization, Hermitian-conjugate conventions and overall
sign are still unverified. Those must be checked before numerical use of this
map. No physical matching status is changed.

## Install and run (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_BC_110101_native_mass_basis_map.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/test_t3_bc_110101_mass_basis_map.py
python -u -X utf8 -m scripts.AuditT3BC110101MassBasisMap --all --run-id bc_110101_map_v1
```

Outputs appear at
`output/audits/t3_bc_110101_mass_basis_map/bc_110101_map_v1/`, as
`B_m1.json`, `C_m1.json` and `summary.json`. Every per-model report carries
source SHA-256s and flags marking **physical matching uncertified**. Choose a
new run ID for reruns; existing results are never overwritten.

## Verification status and next step

- The previous native-quartic mass-basis audit was run locally by the user;
  both models reported `7/7`, rank `2`, and residuals `3.553e-15`, `3.657e-15`.
- This package has passed local Python regression/fixture checks. It still
  requires running against the actual local project before being recorded as
  a completed live-data milestone.
- Next: derive and independently validate the **RGBeta potential prefactors,
  quartic coupling normalization, and conjugation/sign map**, then extend
  source-level proofs to Yukawa and ordinary scalar vertices. Only then attempt
  a native Matchete heavy-mass-basis check and one-loop matching.
- No full sequential EFT matching is performed here.
