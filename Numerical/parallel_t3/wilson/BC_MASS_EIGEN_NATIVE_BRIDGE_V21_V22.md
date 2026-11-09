# Workstream A v21–v22: B/C(−1) native-CG mass-eigenstate transport

## Status and scope

**No physics matching or RGE is performed.** Both stages are read-only and add
new source-conditional JSON diagnostics. The implementation reuses the
existing, source-defined scalar diagonalizer in
`Numerical/experimental/GeneralT3MixedScalarMassBasis.py` instead of
implementing another mass convention.

Source contracts reviewed read-only from repository HEAD
`2bc2e4b516776212645370105cc63495a3c47348`:

- `Numerical/experimental/GeneralT3MixedScalarMassBasis.py`
- `scripts/RunGeneralT3MixedScalarMassBasis.py`
- `scripts/AuditT3BCNativeQuarticMassBasis.py`
- `Numerical/parallel_t3/wilson/complex_scalar_bridge_v19.py` (prior delivered v19)
- `Numerical/parallel_t3/wilson/full_scalar_census_v20.py` (prior delivered v20)
- `Lagrangian/interactions/T3Topology.wl`,
  `Lagrangian/interactions/GeneralScalarOperators.wl`,
  `Lagrangian/model/LagrangianBuilder.wl`.

The implementation uses the *local* source bytes as its acceptance authority;
when any source SHA256 in a saved witness changes, it fails closed.

## What is new relative to existing B/C mass-basis scripts

An existing experimental B/C calculation already diagonalizes the scalar mass
matrix, transports **declaration-level** RGBeta quartics, and separately
rotates the seven historical Stage22 native Matchete extra-quartic tensors.
Those capabilities are not reimplemented here.

**v21** connects four sets of evidence that were previously audited separately:

1. Current native Matchete v18 topology CG: 1 per B/C model.
2. Current native Matchete v18 seven optional extra CG tensors: 7 per model.
3. Existing final candidate B/C scalar mass eigensystems and exact quadratic
   tensor transport, including Hermitian conjugates.
4. v19's explicitly *conditional* six-of-seven baseline RGBeta map and
   v20's 178-slot complete quartic registration census.

It evaluates the exact v18 CG tensors after rotation rather than replacing
them with independently fabricated SU(2) tensors. Each native operator is
checked on random complex field vectors *before and after* transformation,
including its `PlusHc` and an independent complex coefficient. The seventh
optional Matchete direction, which does not belong to the six-coupling
baseline RGBeta special subspace, is explicitly kept and tested nonzero.

**v22** independently rechecks the saved v21 report against its six unchanged
source witnesses, confirms the v20 inventory 150 real + 28 complex and emits
an explicit blocking/remaining-work checklist. It never grants authority to
run physical Wilson evolution.

## Physics and numerical conventions

For B/C(−1) the existing code uses

\[
\Phi=\binom{S_1}{\epsilon S_2^*},\quad
\epsilon=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\quad
\Phi=U\binom{X}{Y},\quad
U^\dagger\begin{pmatrix}m_{S_1}^2&\mu^2\\\mu^{2*}&m_{S_2}^2\end{pmatrix}U
=\operatorname{diag}(m_X^2,m_Y^2).
\]

The inverted field relationship is

\[
S_1=U_{10}X+U_{11}Y,\qquad
S_2=-\epsilon\big(U_{20}X+U_{21}Y\big)^*.
\]

*The actual stored Python implementation uses zero-based row/column indices.*
The `-epsilon` on the `S2` field is essential. For each original native CG
\(T_{abcd}\), with four ordered field legs, the program computes

\[
T'_{ijkl}=T_{abcd}(M_1)_{ai}(M_2)_{bj}(M_3)_{ck}(M_4)_{dl}
\]

for every permitted mass-field choice and sums coincident field monomials.
For arbitrary complex \(c\) and numerical Higgs/eigen-doublet fields,

\[
2\mathrm{Re}[c\,T(S_1,S_2,H)]
=2\mathrm{Re}[c\,T'(X,Y,H)]
\]

is tested at independent points. The quadratic basis is verified against the
existing exact tensor transport and its SU(2) covariance check. The target
relative action tolerance is \(2\times10^{-11}\).

**Caveat on physical signs:** Matchete's `+PlusHc[...]` quartic convention is
being compared on both sides of one *field redefinition*. This does not by
itself identify the minus sign between a potential and its Lagrangian, or the
signed GammaCC→Weyl conversion used in neutrino-mass matching.

## Installation — new files only

All paths are relative to
`C:\thesistemp\T3neutrinomassmodelpipeline`:

1. `Numerical/parallel_t3/wilson/bc_mass_eigen_native_bridge_v21.py`
2. `Numerical/parallel_t3/wilson/bc_mass_bridge_gate_v22.py`
3. `Numerical/parallel_t3/wilson/BC_MASS_EIGEN_NATIVE_BRIDGE_V21_V22.md`
4. `scripts/RunT3ParallelABCMassEigenNativeV21.py`
5. `scripts/RunT3ParallelABCMassBridgeGateV22.py`
6. `scripts/RunT3ParallelABCMassEigenBatchV21V22.ps1`
7. `tests/test_t3_parallel_a_bc_mass_bridge_v21_v22.py`

Install and run (PowerShell):

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$env:USERPROFILE\Downloads\t3_workstream_A_bc_mass_eigen_v21_v22.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
powershell -ExecutionPolicy Bypass -File scripts\RunT3ParallelABCMassEigenBatchV21V22.ps1
if ($LASTEXITCODE -ne 0) { throw "v21/v22 batch failed" }
```

The runner **does not re-run Wolfram**. It uses your already saved evidence:

- `output/parallel_t3/wilson/matchete_topology_native_v18.json`
- `output/parallel_t3/wilson/matchete_bc_extra_native_v18.json`
- `output/parallel_t3/wilson/complex_scalar_bridge_conditional_v19.json`
- `output/parallel_t3/wilson/full_scalar_census_conditional_v20.json`
- `output/general_t3_mass_basis/B_m1_final_mass_basis.json`
- `output/general_t3_mass_basis/C_m1_final_mass_basis.json`

Outputs (created once, independently revalidated on reruns):

- `output/parallel_t3/wilson/bc_native_mass_eigen_complex_conditional_v21.json`
- `output/parallel_t3/wilson/bc_quartic_mass_readiness_conditional_v22.json`

Use the runner's `-VerifyOnly` switch to require existing output files.

## Remaining physical limitations

`MatcheteHeavyMassBasisCheckOnRedefinedLagrangian` remains false. A correct
numerical field transformation and saved native CGs do not constitute native
Matchete's mass-eigenstate Lagrangian, newly registered eigenfields, or the
hard-matched scalar threshold. The numerical native RGBeta triplet component
frame remains independent, and the actual Matchete↔RGBeta physical coupling
map, signed GammaCC→Weyl and physical one-loop Wilson ADM are outstanding.

EFT order stays:

**UV → integrate out ALL fermions → integrate out ALL heavy scalars → SMEFT**.

No fixed-order scalar-threshold leading log may be added again to already
run Wilson coefficients, and no physical Wilson evolution is run by this batch.
