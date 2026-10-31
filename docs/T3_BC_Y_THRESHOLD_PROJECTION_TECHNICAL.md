# B(−1)/C(−1): experimental first scalar-threshold Wilson projection

## Scope and source basis

This package is a **read-only experimental stage** for
`EFT1(H,X,Y) → EFT2(H,X)` in the ordinary T3 B(−1) and C(−1) models. The
existing `scripts/RunT3BCCoupledWilsonSequential.py` was inspected on the
repository's public branch. It uses the saved `continuous_v1` full-scalar
RGBeta trajectory, hard fermion-boundary `C5`, exact saved `LLSS` tensor seed,
selected Stage17 source checks and seven-term coupled Wilson RGE.

The newer `bc_separated_thresholds_v1` study described in the user's handover
uses a **different candidate scalar trajectory**. Its local-only script and
JSON report were unavailable to inspect here. Consequently this package **does
not import the standalone threshold numbers**. It recomputes an internally
consistent running Y threshold directly on the same continuous-scalar state
that drives the Wilson ODE. It does not assert agreement with the standalone
numbers (`B: 1187.55724 GeV`, `C: 1185.14817 GeV`).

No production source file or previous result is changed. The package records
SHA256 hashes of every file it reads, together with its script, numerical
engine and mass-basis module. Endpoints are cross-checked against existing
saved continuous-study checkpoints. The Wilson engine's preexisting grouped
fermion hard boundary and direct/Stage17 source parity checks are retained.

## Exact field transformation

All complex doublet components use

\[
 z_a=\frac{r_{2a}+i r_{2a+1}}{\sqrt2}.
\]

The original and eigenstate fields obey

\[
\begin{pmatrix}S_1\\ \epsilon S_2^*\end{pmatrix}
=U\begin{pmatrix}X\\Y\end{pmatrix},\qquad
\epsilon=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\qquad
U^\dagger A U=\operatorname{diag}(m_X^2,m_Y^2),
\]

where

\[
A=\begin{pmatrix}m_{S_1}^2&\mu_2\\ \mu_2^*&m_{S_2}^2\end{pmatrix}.
\]

In particular

\[
S_{1,a}=\sum_{k=X,Y}U_{1k}\chi_{k,a},\qquad
S_{2,a}=-\sum_{k,b}\epsilon_{ab}U_{2k}^*\chi_{k,b}^*.
\]

If `old = T @ new` is the corresponding **real** field map, the real
2×2 blocks for a complex coefficient \(z\) are

\[
R(z)=\begin{pmatrix}\operatorname{Re}z&-\operatorname{Im}z\\
\operatorname{Im}z&\operatorname{Re}z\end{pmatrix},\quad
R_*(z)=\begin{pmatrix}\operatorname{Re}z&\operatorname{Im}z\\
\operatorname{Im}z&-\operatorname{Re}z\end{pmatrix}.
\]

The nontrivial 2×2 real blocks of \(T\) are

\[
T_{S_{1,a},\chi_{k,b}}=\delta_{ab}R(U_{1k}),\qquad
T_{S_{2,a},\chi_{k,b}}=-\epsilon_{ab}R_*(U_{2k}^*).
\]

Each listed block acts on the `(Re, Im)` pair. The Higgs block is identity.
The actual implementation constructs the full real \(12\times12\) matrix
and checks \(T^TT=TT^T=I\), independently reproduces the complex field
identity and verifies intertwinement under random SU(2)×U(1) transformations.
The hypercharges in the stated \(Q=T_3+Y\) convention are
`Y(H)=Y(S2)=+1/2` and `Y(S1)=Y(X)=Y(Y)=-1/2`.

## Wilson coefficient formula and finite matching limitation

The existing Wilson array uses six independent physical lepton flavour pairs
`(00,01,02,11,12,22)` and shape

\[
 C^{\mathrm{old}}_{ijab,f}\in\mathbb C^{3\times3\times12\times12\times6}
\]

with fermion gauge-index basis `L1,L2,eC`. Because the scalars are
**real-valued component fields**, the transformation is a transpose, not a
Hermitian conjugation:

\[
 C^{XY}_{ijAB,f}=\sum_{a,b=1}^{12}
 T_{aA}T_{bB}\,C^{\mathrm{old}}_{ijab,f}.
\]

Here the fermion basis and six-flavour channel are **not** rotated. The output
at the high-scalar running crossing is

\[
 C^{\mathrm{EFT2,tree}}_{ijAB,f}
 = C^{XY}_{ijAB,f}\big|_{A,B\in\{H,X\}}
 \in\mathbb C^{3\times3\times8\times8\times6}.
\]

The projected components include `HH`, `HX` (with the symmetric `XH`
partner), and `XX` blocks, with lepton flavour and fermion gauge indices
unchanged. The output records their norms and occupied-component counts;
no invented operator coefficient is inserted. For a subsequent **EFT2
Wilson ODE**, the scalar-dependent tensors in the existing seven-term RGE
must be constructed again using *only* the real components of `H,X`, the
remaining EFT2 quartics and gauge generators, and the matching of the other
renormalizable couplings. Reusing the 12-scalar kernels unchanged would be
incorrect. The existing one-generation SM anomalous-dimension limitation
persists. Further operators or finite thresholds, if required by the EFT2
matching calculation, have not been derived here.

The full Wilson threshold has the structure

\[
 C_{\mathrm{EFT2}}(\mu_Y^-)
 = C_{\mathrm{EFT2,tree}}(\mu_Y^+)
   +\Delta C_Y^{\mathrm{finite,1\!\! -\!\! loop}} .
\]

**The second term is unknown, not zero.** This package does not calculate
that term, does not run EFT2, does not locate the EFT2-corrected X threshold,
and does not calculate a physically matched neutrino mass spectrum.

The `HH` Weinberg coefficient is recorded **only as a diagnostic**. Since the
Higgs is unchanged under scalar diagonalisation, its value and complete
30-component SU(2) Weinberg-template residual must be invariant under the
rotation. This is not a claim that final `C5` at low energy has been matched.

## Numerical crossing and archived output

The continuous RGBeta source is rerun to obtain `mS1Sq(t)`, `mS2Sq(t)`,
`mu2(t)`, all scalar quartics and gauge/Yukawa backgrounds. Within an
adjacent pair of source nodes, the quadratic entries are interpolated
**linearly in** \(t=\ln\mu\). The root

\[
\lambda_{\rm max}(A(t_Y))=\exp(2t_Y)
\]

is located with Brent's method and required to be unique on the saved EFT1
range. The complete original-basis Wilson ODE is integrated only down to
\(\mu_Y=\exp t_Y\), with its last kernel interpolated using the Wilson
engine's existing logarithmic-grid kernel interpolator.

For each model, the new run produces:

- `<tag>_y_threshold_wilson.json`: provenance, numerical crossing, scalar
  unitarity/covariance, field/tensor identity, pair-symmetry checks,
  historical source parity, survivor norms and **false** physical matching
  certification flags.
- `<tag>_y_threshold_wilson.npz`: three complete complex tensor arrays
  `C_EFT1_at_Y_original_real`, `C_EFT1_at_Y_XY_real`,
  `C_EFT2_tree_at_Y_HX_real`, and the real mixing matrix, complex U,
  threshold and ODE scale nodes.

The NPZ is the handoff to subsequent EFT2 work; do not replace it with the
six projected Weinberg coefficients.

## Verification status

Completed in this package's standalone development environment:

- Python syntax compilation of the new experimental module and CLI runner.
- 11 NumPy-only unittest regressions, including negative controls for
  the epsilon sign, complex mixing phase, old-field ordering, tensor index
  ordering and wrong SU(2) scalar action.

**Not performed:** Windows-local RGBeta re-integration, continuous-study
source-parity validation with saved data, Matchete, and a numerical B/C
end-to-end \(\mu_Y\) run. Those require the up-to-date local repository.

## Installation and execution

The ZIP contains complete new files at these **project-relative paths**:

- `Numerical/experimental/T3BCScalarYThresholdProjection.py`
- `scripts/RunT3BCYThresholdProjection.py`
- `tests/test_t3_bc_y_threshold_projection.py`
- `docs/T3_BC_Y_THRESHOLD_PROJECTION_TECHNICAL.md`
- `docs/T3_BC_Y_THRESHOLD_TREE_PROJECTION_MILESTONE.md`

From PowerShell at `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
python -m unittest discover -s tests -p "test_t3_bc_y_threshold_projection.py" -v
python -u -X utf8 -m scripts.RunT3BCYThresholdProjection --all --run-id continuous_v1 --output-id bc_y_threshold_projection_v1
```

The default output directory is
`output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/bc_y_threshold_projection_v1/`.
The prerequisites are the existing `direct_wilson_sequential_v1` and
`bc_extra_wilson_sequential_v1` reports under the same `continuous_v1`
run folder, each with the same `--samples` default of 65. If the user-local
names or files differ, use the CLI's `--direct-output-id`,
`--selected-output-id`, `--source-root` and `--run-id` arguments rather than
fabricating data or silently using the other threshold trajectory.

**Resolution log:** The local authoritative
`docs/T3_GENERAL_SCALAR_MATCHING_RESOLUTION_LOG.md` is not included as a
replacement because its recent changes may exist only on Windows. Append
`docs/T3_BC_Y_THRESHOLD_TREE_PROJECTION_MILESTONE.md` once to that file after
review rather than overwriting local milestones.
