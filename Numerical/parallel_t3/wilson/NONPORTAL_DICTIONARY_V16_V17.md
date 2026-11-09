# Workstream A v16–v17 — all 76 native Matchete non-portal real quartics

## Verified source contracts (public main, read-only, 2026-11-09)

- `Lagrangian/interactions/ScalarPotential.wl`: `BuildScalarSelfInvariantCandidates`, `BuildCrossScalarInvariantCandidates`, `PhysicalSU2InvariantBasis`; Matchete real coefficients of scalar self have a **1/2** overall factor, scalar cross has **1**. Native CGs include identical-boson symmetrization.
- `scripts/RunT3ParallelAExportPortalMatchete.wl`: provides the exact native source status gates, CG accessor, and B(-1)/C(-1) *symbolic-only* exception. The new v16 exporter uses the same explicit status convention and never calls `Match[]`.
- `RGE/running/rgbeta/T3RGBetaModel.wl`: `lambdaS1`, `lambdaS2` are norm-type self invariants each carrying 1/2; triplet self-adjoint `lambdaS1Adj`, `lambdaS2Adj` also carry 1/2; `lambda12`, `lambda12Adj`, and the E-class `lambda12Cross` use their source factors **1**. No changes to this file.
- RGBeta v14 native PackageScope `Lam`/projector has been run by the user: **150/150** real diagonals and **110/110** Higgs-portal response entries. v15 complex channels are **exploratory** only.

## v16 native evidence and conditional action audit

Native Matchete export accesses `GetCGTensor` for all self- and cross-scalar CGs in all 16 models, with a source SHA256 witness and explicit `LInt`/allowed-interaction checks. It writes:

`output/parallel_t3/wilson/matchete_nonportal_native_v16.json`

The Python module checks current source hashes for this export, the v5 native Higgs portals, and the independently executed native v14 RGBeta projector witness. It reconstructs the scalar representation matrices **from the Matchete portals**, confirms their SU(2) Lie brackets, Casimir and trace normalization, and fits native Matchete CGs to the source-defined RGBeta invariant basis.

For a triplet, the source-conditioned bilinear adjoint is recovered from the Higgs portal:

\[
K_{ijab}=\sum_A(t_H^A)_{ij}(T_S^A)_{ab},\qquad
(T_S^A)_{ab}=2\sum_{ij}(t_H^A)^*_{ij}K_{ijab}.
\]

For non-portal scalar tensors, the target invariants are the identity, \(\sum_A T^A_1\otimes T^A_2\), and, in the joint triplet–triplet sector, the transported crossed invariant. The Matchete triplet may use a **spherical/weight basis**, not the real Cartesian basis implicit in RGBeta's `del[S2]`. The implementation derives the invariant symmetric bilinear form `Q` independently from the portal-reconstructed generators by solving `T_a^T Q + Q T_a = 0`, normalizes `Q† Q = I`, and constructs `conj(Q)_{ac} Q_{bd} + δ_{ad} δ_{bc}`. In a Cartesian frame `Q = I`, so this reduces to RGBeta's `δ_{ac}δ_{bd} + δ_{ad}δ_{bc}`. This avoids assuming Cartesian axes in native Matchete data. Scalar-self tensors are symmetrized on two barred and two unbarred identical-boson legs before the fit. All local Matchete-to-RGBeta matrices must be real and full rank. Random independent complex scalar actions test the reconstructed polynomial, not just the tensor fit.

Counts across 16 models:

| Family | Self S1 | Self S2 | Cross S1/S2 | Total |
|---|---:|---:|---:|---:|
| A (3 points) | 3 | 6 | 3 | 12 |
| B (3 points) | 3 | 3 | 6 | 12 |
| C (3 points) | 3 | 3 | 6 | 12 |
| D (3 points) | 6 | 3 | 3 | 12 |
| E (4 points) | 8 | 8 | 12 | 28 |
| **Total** | **23** | **23** | **30** | **76** |

The numbers are *expected source registrations*, not claimed successful native Matchete evaluations until your Windows export completes. Output:

`output/parallel_t3/wilson/nonportal_dictionary_conditional_v16.json`

## v17 completeness cross-check

A fresh audit reruns the v16 mathematical/source checks **and** the v14 source-guarded audit (which in turn rechecks v9 including native v5–v8 source evidence). The 150 real slots divide as:

- 16 SM Higgs norm terms;
- 58 Higgs–scalar portal directions from v9;
- 76 non-portal scalar self/cross directions from v16.

Native RGBeta's 150 real projector diagonal evaluations and 110 portal response entries are independently checked. This is the real-sector **registration and conditional action dictionary census**, not a unique physical scalar-basis map. Output:

`output/parallel_t3/wilson/full_scalar_census_conditional_v17.json`

## Limitations (cannot promote to physical)

1. The S2 frame is *reconstructed from Matchete Higgs portals*. RGBeta's numerical S2 frame is not separately measured. A basis-equivalent tensor/action relation is verified only conditionally.
2. B(-1)/C(-1) are **not** in the heavy scalar mass eigenbasis. Native CG extraction is allowed; mass matching is not.
3. The 28 RGBeta complex quartic channels still lack a complete **native Matchete complex action** certificate; the pre-existing 50 extra-quartic dictionary and the v11/v12 constructed-action checks are separate conditional evidence.
4. Signed GammaCC→Weyl and the complete one-loop Wilson ADM are unverified. No physical RGE or threshold evolution may be authorized by v16/v17.

The EFT order remains **UV → integrate out ALL fermions → integrate out ALL scalars → SMEFT**.

## Installation and commands (PowerShell)

Download the ZIP and expand it into `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$env:USERPROFILE\Downloads\t3_workstream_A_nonportal_v16_v17.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
powershell -ExecutionPolicy Bypass -File scripts\RunT3ParallelAScalarDictionaryV16V17.ps1
```

The runner performs both test modules, creates the native v16 Matchete CG JSON *only when absent*, verifies it through the conditional v16 audit, and revalidates existing v16/v17 outputs instead of overwriting them. Failure preserves existing evidence and exits nonzero. No source Lagrangian, RGBeta installation, Matchete installation, B/C workstream files, or GitHub repository is modified.

**Expected if the native 76-slot tensor fit passes:** 76/76 Matchete nonportal CG directions, 48/48 self/cross sectors, and 150/150 conditional real-quartic census. These remain conditional results.
