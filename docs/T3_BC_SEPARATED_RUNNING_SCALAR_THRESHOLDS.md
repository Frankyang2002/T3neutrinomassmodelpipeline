# B/C(−1): running Y/X threshold preparation (candidate EFT1)

## Scope

This is the next **numerical preparation stage** toward the physical
`F → Y → X` sequence. It locates the running scalar mass thresholds, retains
all mixing phases, produces scalar-flavour projectors and the corresponding
Matchete-native Yukawa weights. It does **not** perform Y/X Wilson matching,
finite one-loop threshold corrections, or EFT2 running after removing Y.

It **uses the same `output/general_t3_quadratic_eft1` benchmark as the
successful saved one-scale C5 replay**. Do not substitute the separate
`continuous_v1` run: that run has different quadratic mass inputs.

## Equations and conventions

The two hypercharge-−1/2 doublets form

\[
\Phi=\begin{pmatrix}S_1\\\varepsilon S_2^*\end{pmatrix}
=U(\mu)\begin{pmatrix}X\\Y\end{pmatrix},\qquad
\mathcal{M}^2(\mu)=\begin{pmatrix}
 m_{S_1}^2(\mu)&\mu_2(\mu)\\
 \mu_2^*(\mu)&m_{S_2}^2(\mu)
\end{pmatrix}.
\]

The code uses the project `MixedDoubletMassBasis` class to diagonalise this
Hermitian matrix. It then solves separately

\[
\mu_Y=\sqrt{\lambda_+(\mu_Y)},\qquad
\mu_{X,\mathrm{un-decoupled}}=\sqrt{\lambda_-(\mu_{X,\mathrm{un-decoupled}})},
\]

using Brent root finding in \(t=\ln\mu\), with an actual DOP853 integration of
the saved candidate EFT1 beta functions. The initial and final state, all
additional quartic beta keys, exact source beta hash, and saved one-scale
mass-basis inputs must agree; otherwise the run fails.

The native Matchete scalar Yukawa factors are

\[
Y_a^{(1)}=U_{1a}y_1,\qquad
Y_a^{(2)}=-\alpha_2 U_{2a}y_2,
\]

with \(\alpha_2=-1\) for B(−1) and \(\alpha_2=+1\) for C(−1), as
verified by the previous native Matchete certificates. The source vertices
have distinct fermion conjugations and must **not** be merged into a single
Yukawa operator. The rank-one mixing-space threshold projector is
\(P_a=U_{:a} U_{:a}^{\dagger}\) in the \(\Phi\) basis.

## Verification gates

1. Saved `RunGeneralT3QuadraticEFT1` report and exact hashed beta source.
2. Native Matchete `HeavyMassBasis` and distinct Yukawa CG certificates.
3. Identical quadratic-EFT1 benchmark to the one-scale C5 report (SHA and
   independent endpoint mass matrix parity).
4. Recomputed EFT1 DOP853 endpoint matches its saved mass, mixing and extra
   quartic values.
5. Two positive, ordered mass-crossing roots in
   `1000 GeV < X < Y < 100000 GeV`; projectors are orthogonal.

## Files

- `scripts/RunT3BCSeparatedScalarThresholds.py` — new read-only runner.
- `tests/test_t3_bc_separated_scalar_thresholds.py` — numerical regression tests.
- `docs/T3_BC_SEPARATED_RUNNING_SCALAR_THRESHOLDS.md` — this note.
- `docs/T3_GENERAL_SCALAR_MATCHING_RESOLUTION_LOG.md` — progress log.

### PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_BC_separated_scalar_thresholds_v1.zip' -DestinationPath . -Force
python -X utf8 -m pytest -q tests/test_t3_bc_separated_scalar_thresholds.py
python -u -X utf8 -m scripts.RunT3BCSeparatedScalarThresholds --all --run-id bc_separated_thresholds_v1
```

Reports appear in:
`output/audits/t3_bc_separated_scalar_thresholds/bc_separated_thresholds_v1/`.
Each model has a complete per-root rotation/projector/vertex-weight record.
The output does **not** add or replace saved candidate neutrino masses.

## What remains for full sequential matching

- Start from the existing heavy-F hard Wilson boundary. The new script does
  not recompute or endorse the hard coefficient.
- Evolve the *matched* full LLSS Wilson tensor down to \(\mu_Y\), with
  compatible mass-eigenstate gauge-component conventions.
- Derive the actual Y-removal operator map and finite matching coefficients.
- Implement the **Y-absent, X-active EFT2** beta functions and run the
  remaining Wilson tensor to the now shifted \(\mu_X\).
- Match out X with the finite scalar correction; then run SM+Weinberg C5.

The reported X scale is **explicitly only the un-decoupled EFT1 estimate**:
continuing the two-scalar beta system below \(\mu_Y\) is not the actual
post-threshold EFT. Do not describe the output as a complete sequential
calculation or physical neutrino prediction.
