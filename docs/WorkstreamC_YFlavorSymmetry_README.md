# Workstream C v8 — Native Y/LLHH flavour symmetry and factorization

**Scope:** exact native *neutral SU(2) component* with `Proj[1]` (PR), its `Proj[-1]` (PL) Hermitian conjugate, and formal finite/pole Yukawa flavour structure. This is **not** a physical Weinberg coefficient, a `GammaCC` → two-component Weyl normalization, a full EFT renormalization, or an oscillation-parameter prediction.

## Verified input provenance

The source files are existing repository outputs (read-only inspected before writing this utility):

- `output/audits/t3_bc_y_renormalization_ledger_v1/B_m1_y_llhh_constant_neutral_flavour.json`
- `output/audits/t3_bc_y_renormalization_ledger_v1/B_m1_y_llhh_pole_neutral_flavour.json`
- `output/audits/t3_bc_y_renormalization_ledger_v1/C_m1_y_llhh_constant_neutral_flavour.json`
- `output/audits/t3_bc_y_renormalization_ledger_v1/C_m1_y_llhh_pole_neutral_flavour.json`

The input LEDGER files express the **ordered** neutral coefficient in `hbar/MF` units, with Higgs SU(2) component 2 and lepton component 1 and `epsilon12=+1`. The source-level PR/PL counts are 16/16 for formal finite and 25/25 for pole. Each PR monomial contains a quartic, four mixing spurions, and a *pair of unbarred* Yukawas `y1`/`y2`; the PL counterpart contains exactly conjugated couplings.

### New exact result: flavour structure

Define the 3-vectors of Yukawa couplings `y1`, `y2`, and the `3 x 2` matrix `Y=[y1 y2]`. The native PR **ordered** flavour kernel has the form

\[
 A_{ij}^{\rm ord}=a_{11}(y_1)_i(y_1)_j+a_{12}(y_1)_i(y_2)_j+a_{22}(y_2)_i(y_2)_j,
\]

where `a11`, `a12`, `a22` are *exact symbolic expressions in the three Higgs quartics and mass-eigenstate mixing spurions*. The ledger has **no separately ordered** `y2_i y1_j` terms; do not invent one in the source.

As a **pure algebraic flavour projection**, not a Wilson normalization,

\[
 A^{\rm sym}=\tfrac12(A^{\rm ord}+A^{{\rm ord}\,T})
 =Y\begin{pmatrix}a_{11}&a_{12}/2\\a_{12}/2&a_{22}\end{pmatrix}Y^T.
\]

This is symmetric under **transpose, not Hermitian conjugate**. For two independent vectors, its rank is at most two. Explicitly `n=y1 cross y2` obeys `A_sym n=0` (over complex flavour space, using transpose). The same rank bound holds for this isolated native pole and formal finite scalar-bubble sector; it does **not automatically certify a massless neutrino in the full physical T3 calculation**, because all matched pieces, RG effects and normalizations are not yet assembled.

**Important:** Do not set `C5 = A_sym` or `C5 = 2*A_sym` merely based on this construction. Whether the physical convention `L5 = (1/2) C5 O + h.c.` requires a separate overall factor and sign cannot be fixed by PR/PL conjugacy alone. Do not sum PR and PL; they represent Hermitian-conjugate chiral sectors.

### New exact result: B/C relative group signs

For *identical formal* scalar quartic and mixing-spurion names, the source ledgers satisfy **in both pole and formal finite sectors**:

\[
 a_{11}^{C}=a_{11}^{B}/\sqrt3,\quad
 a_{22}^{C}=a_{22}^{B}/\sqrt3,\quad
 a_{12}^{C}=-a_{12}^{B}/\sqrt3.
\]

This relative sign is a verified **native expression** property, not an equality of numerical physical benchmark Yukawas or quartic couplings. The actual B and C scalar mixing and couplings can differ; the script does not silently equate their benchmark values. The cross-Yukawa flip is important when symmetrizing the kernel and testing neutrino textures.

## Files

Only one complete new Python script and this README are added, preserving project paths:

- `scripts/AnalyzeT3BCYFlavorSymmetry.experimental.py`
- `docs/WorkstreamC_YFlavorSymmetry_README.md`

All inputs are existing files. The script uses Python `sympy` (already installed in your project venv), exact rational/source-monomial arithmetic, strict schema checks, exact PL/PR Hermitian conjugacy, explicit matrix factorization, and a self-test for rank and transpose. It writes the actual symbolic B/C coefficients and their SHA-256 input hashes to JSON. It deliberately does not perform Mathematica matching, modify the source Lagrangian or apply any counterterm.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "$HOME\Downloads\T3_WorkstreamC_YFlavorSymmetry_v8.zip" `
  -DestinationPath . -Force

python scripts\AnalyzeT3BCYFlavorSymmetry.experimental.py --self-test

python scripts\AnalyzeT3BCYFlavorSymmetry.experimental.py `
  --output output\audits\t3_bc_y_flavour_symmetry_v1\summary.json
```

Expected: both benchmarks report exact `PR/PL` conjugacy and rank `<=2` factorization, with C-versus-B mixed-Yukawa sign reversal. The output is **not** a mass matrix.

## Next genuine physics step

Verify the `GammaCC`/projector bilinear-to-Weyl normalization in the installed native Matchete version (ideally by independently expanding one component in a fixed charge-conjugation convention or comparing a known Majorana operator's Feynman rule). Then establish the full EFT1/EFT2 one-loop operator and field counterterm basis, and only then use the resulting physical `C5` in sequential Y→X running and mass extraction. The v7 result indicates that the neglected finite-mass corrections from `r=MX^2/MY^2 ~ 0.76` are not small; do not incorporate them by multiplying the physical Wilson coefficient by the full-mass bubble without proper higher-dimensional matching.
