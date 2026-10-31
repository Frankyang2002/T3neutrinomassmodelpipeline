# B(−1)/C(−1): native Matchete extra-quartic mass-basis audit

## Purpose and status

This is the first source-tensor verification step toward resolving the B/C(−1)
mixed-scalar mass-basis obstruction. It does **not** resolve the complete
Matchete/RGBeta matching problem. The new audit uses the **actual saved native
Matchete tensors**, rather than the arbitrary-source-tensor identity tests
already present in `RunGeneralT3MixedScalarMassBasis`.

The source inputs must exist:

- `output/general_t3_stage22/matchete_tensors.json` (native Matchete quartics)
- `output/general_t3_mass_basis/B_m1_final_mass_basis.json`
- `output/general_t3_mass_basis/C_m1_final_mass_basis.json`

These are read-only. No previously fitted parameters, RG trajectories,
Lagrangian definitions, or matching outputs are replaced.

## What is checked

For B/C(−1), define the hypercharge \(-1/2\) doublet pair and mass eigenfields:

\[
\Phi=\begin{pmatrix}S_1\\\varepsilon S_2^*\end{pmatrix}
     =U\begin{pmatrix}X\\Y\end{pmatrix},\qquad
\varepsilon=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\quad
U^\dagger\mathcal M_\Phi^2 U=\operatorname{diag}(m_X^2,m_Y^2).
\]

Each native four-field component tensor \(T_{abcd}\) is transformed using

\[
T'_{ijkl}=T_{abcd}(M_1)_{ai}(M_2)_{bj}(M_3)_{ck}(M_4)_{dl},
\]

where each \(M_k\) is taken from the existing
`GeneralT3MixedScalarMassBasis.field_map`, including the essential
\(-\varepsilon\) action for \(S_2\) and \(S_2^\dagger\).

The program checks:

1. All **seven** B/C native extra-quartic invariant tensors are present,
   have the expected ordered fields and complete \(2^4\) components.
2. Both original and rotated quartic polynomials agree at random complex
   field configurations, including a complex coupling times the Hermitian
   conjugate; tolerance \(3\times10^{-12}\).
3. The two independent Matchete invariants with field counts `110101`
   retain **rank 2** for a common rotated field signature. Neither tensor
   is replaced by the other or silently discarded.
4. The saved scalar mass-matrix rotation agrees with the existing
   `MixedDoubletMassBasis` diagonalization.
5. Input SHA-256 hashes and the full rotated native component tensors are
   saved for reproducibility.

Only the seven **extra** Matchete quartics are tested by this package.
The original scalar quartics and Yukawa Clebsch–Gordan tensors have not been
rotated from native sources in this audit.

## Install and run (PowerShell)

Place the ZIP in a known location and extract from the local repository root:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_BC_native_quartic_mass_basis_audit.zip' -DestinationPath . -Force
python -X utf8 -m pytest -q tests/test_t3_bc_native_quartic_mass_basis.py
python -u -X utf8 -m scripts.AuditT3BCNativeQuarticMassBasis --all --run-id bc_native_vertices_v1
```

If the source mixed-basis reports do not exist locally, generate them first:

```powershell
python -u -X utf8 -m scripts.RunGeneralT3MixedScalarMassBasis --all
```

Use a new `--run-id` for reruns; existing audit output is not overwritten.

New outputs:

`output/audits/t3_bc_native_quartic_mass_basis/bc_native_vertices_v1/`

- `B_m1_native_quartic_mass_basis.json`
- `C_m1_native_quartic_mass_basis.json`
- `summary.json`

## What is resolved, and what is still open

**Resolved if these audits pass:** The seven saved native Matchete extra-quartic
component tensors can be transported through the numerical B/C mixed-doublet
mass-basis map without losing a tensor or a second independent contraction.
This improves upon the previous *abstract* tensor-transport check.

**Still open, explicitly not certified by this audit:**

- Native *baseline* quartic, Yukawa and kinetic/gauge-field component vertices
  must be consistently transformed and checked.
- RGBeta native invariant tensors, normalization and coupling-sign/basis maps
  must be independently compared with the Matchete tensors. In particular,
  the two `110101` contractions must be mapped as a **2×2** coupling-basis
  transformation, not a guessed scalar sign.
- A newly assembled Matchete Lagrangian **in the mass basis** must pass
  `CheckLagrangian` and its `HeavyMassBasis` diagnostic.
- Independent mass-eigenstate one-loop matching and then physically consistent
  fermion/scalar sequential EFT matching remain to be certified.

The existing B/C one-scale numerical neutrino-mass results remain **conditional**.
Do not change the `physical_*_certified=False` flags based on this audit alone.

## Change record

- Added `scripts/AuditT3BCNativeQuarticMassBasis.py`: read-only saved-native
  tensor rotation with full component outputs and fail-closed evidence.
- Added `tests/test_t3_bc_native_quartic_mass_basis.py`: native tensor,
  complex-bilinear, two-contraction and malformed-source regression tests.
- Added this `docs/T3_BC_NATIVE_QUARTIC_MASS_BASIS_AUDIT.md` record.
- No existing physics or production files modified.
- All assertions require local verification; this file does not assert that
  live Matchete or the user's local benchmark was run by the package author.
