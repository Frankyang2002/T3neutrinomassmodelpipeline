# Workstream A v6 — SU(2) portal adjoint orientation

## Scope

New-only A-owned files, no changes to current scalar potential, RGBeta model,
Wilson solver, matching, historical reports, or GitHub. v5's user-confirmed
58/58 Matchete CG extraction, 58/58 RGBeta portal registrations and 14 corrected
anti-Hermitian `I` phases are the input evidence.

Two levels of evidence must **not** be conflated:

1. **Native RGBeta SU(2) fundamental Fierz identity** (directly testable in
   Mathematica):

   \[
   \sum_A (t_f^A)_{ij}(t_f^A)_{kl}
   =\tfrac12\left(\delta_{il}\delta_{kj}
                 -\tfrac12\delta_{ij}\delta_{kl}\right),
   \quad T(f)=1/2.
   \]

   For native Matchete portal doublets sharing the fundamental index axes, this
   fixes the adjoint **relative sign**, not merely its norm. For the registered
   B/C doublet CG fixture, ordered as `(lambdaH{n}Inv1,lambdaH{n}Inv2)`, the
   source-action coefficients are

   \[
   \binom{\lambda_{Hn}}{\lambda_{Hn\mathrm{Adj}}}
   =\begin{pmatrix}1/2&\sqrt3/2\\2&-2/\sqrt3\end{pmatrix}
    \binom{\lambda_{Hn\mathrm{Inv1}}}{\lambda_{Hn\mathrm{Inv2}}}.
   \]

   The v6 Python audit checks this matrix from **actual local v5 exports**, not
   from coupling names, and independently holds out random complex potentials.
   The historical five-class fixture is not a substitute for the live run.

2. **Triplet Pauli-frame Lie-algebra check:** RGBeta's `SU2L[S2]` generators
   remain abstract (`tGen` has no numeric component frame). Nonetheless the
   *Matchete registered Cartesian* portal CG, after `I`-phase Hermitian
   projection, defines a traceless invariant `K`. With Higgs
   \(t_H^a=\sigma^a/2\), reconstruct

   \[
   (G_R^a)_{uv}=2\sum_{ij}(t_H^a)_{ji}K_{ijuv},\quad
   \lVert K\rVert_F^2=3.
   \]

   Only one of the candidate signs can satisfy
   \([G_R^1,G_R^2]=iG_R^3\). Require also Hermiticity, zero trace,
   \(\mathrm{tr}(G_R^aG_R^b)=2\delta^{ab}\), Casimir \(\sum_aG_R^aG_R^a=2I_3\),
   and covariance of `K` under the joint generators. The **v5 pivot-defined
   orientation** for the registered five-class fixture has the wrong Lie
   bracket sign; v6 selects its negative. This is a **Pauli/Matchete-frame
   statement**, NOT a directly extracted RGBeta S2 generator matrix. RGBeta's
   source only specifies the abstract S2 symmetric-square contraction and
   Casimir/trace normalizations. An independent native cross-package S2
   component intertwiner remains missing.

## What remains blocked

- Numeric native RGBeta `SU2L[S2]` generator components in the same Matchete
  Cartesian coordinates and the independent scalar-coordinate intertwiner.
- Full baseline self/cross quartic and scalar-mass dictionary; B(-1)/C(-1)
  charged-scalar mixed mass diagonalization.
- Yukawa and signed GammaCC -> Weyl normalization; validated Wilson ADM and
  scalar thresholds. The mandatory EFT ordering is:
  **UV → integrate out all fermions → integrate out all scalars → SMEFT.**

**Do not load a v6 conditional triplet coefficient directly into the physical
Wilson background.** The result is a component/convention check only.

## Installation and run (PowerShell from repository root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "C:\path\to\t3_workstream_A_portal_orientation_v6.zip"
Expand-Archive -Path $zip -DestinationPath . -Force

.\.venv\Scripts\python.exe -m unittest discover -s tests -p "test_t3_parallel_a*.py" -v

wolframscript -file scripts\RunT3ParallelAExportNativeFierz.wl `
  output\parallel_t3\wilson\rgbeta_fundamental_fierz_native_v6.json

.\.venv\Scripts\python.exe -m scripts.RunT3ParallelAOrientPortals `
  --matchete output\parallel_t3\wilson\matchete_portals_native_v5.json `
  --rgbeta output\parallel_t3\wilson\rgbeta_portals_native_v5.json `
  --v5 output\parallel_t3\wilson\portal_dictionary_conditional_v5.json `
  --fierz output\parallel_t3\wilson\rgbeta_fundamental_fierz_native_v6.json `
  --project-root . `
  --output output\parallel_t3\wilson\portal_orientation_conditional_v6.json
```

Expected after live native witness succeeds:

```
12/12 doublet Fierz, 14/14 triplet Pauli-Lie, 6/6 singlets
DoubletNativeRGBetaFierzFixed_TripletPauliLieSignFixedConditional_NativeS2IntertwinerUnverified
```

All inputs are immutable; each output is **new and refused if already exists**.
The Python audit source-checks the v5 native exporters and the new live Fierz
witness and re-runs the v5 tensor audit before assigning signs.

## Files supplied

- `Numerical/parallel_t3/wilson/portal_orientation_v6.py`
- `Numerical/parallel_t3/wilson/PORTAL_ORIENTATION_V6.md`
- `scripts/RunT3ParallelAExportNativeFierz.wl`
- `scripts/RunT3ParallelAOrientPortals.py`
- `tests/test_t3_parallel_a_portal_orientation_v6.py`
