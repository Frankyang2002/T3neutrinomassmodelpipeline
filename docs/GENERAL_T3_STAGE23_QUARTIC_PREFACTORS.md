# Stage 23: native quartic prefactors and conditional signed coupling maps

## Scope

This is a **read-only diagnostic** on the Stage 22 post-fix native evidence. It does **not** modify `T3RGBetaModel.wl`, beta functions, numerical inputs, or any Stage 22 certification flags. The result is a *conditional signed map*, not a certified physical coupling transformation. The calculation fails closed if the `B_-1`/`C_-1` corrected `110101` invariant is absent or if the Stage 22 exports are stale/incomplete.

## Source-verifiable RGBeta factors

The RGBeta `Package/FieldsAndCouplings.m` implementation (upstream `aethomsen/RGBeta`) uses, in `AddQuartic`:

```wolfram
normalization = 24;
Do[
  If[!$scalars[phi, SelfConjugate] || Head @ phi === Bar,
    normalization *= 1/Sqrt[2]
  ],
  {phi, {phi1, phi2, phi3, phi4}}
];
lam = Evaluate[normalization Matrix[coupling] @ ##] &;
```

Consequently, for four complex scalar slots, the registered quartic tensor coefficient is

\[
n_{\mathrm{RGBeta}} = 4!\left(\frac{1}{\sqrt2}\right)^4 = 6.
\]

The native `rgbeta_native_declarations.json` and `rgbeta_native_invariants.json` files are checked for the **actual registered** `6*Matrix[lambda][##1]&` on **all 50 slots**, including the reuse cases. `SelfConjugate` must be `False` and the two independent exports must agree on their invariant functions and permutation lists.

In the same upstream source, `CouplingTensorMap` calls `AvgPermutations`, assigning each of the `N` `UniqueArrangements` terms weight `1/N`. With the standard real-field quartic contraction `\Lambda_{abcd}\phi_a\phi_b\phi_c\phi_d/4!`, their total coefficient per unsymmetrised complex-field orientation is

\[
\frac1{4!}\sum_{p=1}^N\frac6N=\frac14
=\left(\frac1{\sqrt2}\right)^4.
\]

The matching `1/4` is a consistency check of the **complex-to-real component normalisation**. In particular, **the visible `6*Matrix[...]` is not an extra factor of six in the physical scalar potential**. This check does not independently identify RGBeta's live real-scalar basis with Matchete's basis and does not establish the physical relative sign.

The Stage 22 Matchete tensor exporter marks all 50 quartic contributions with `potential_prefactor="1"`, and the Matchete `Lagrangian/interactions/ScalarPotential.wl` constructs the scalar interactions as Lagrangian expressions, with `PlusHc` for non-Hermitian terms. This annotation is source evidence, not a proof that all component conventions have been globally identified.

## Exact conditional coupling map

For each field-content block, the Stage 22 aligned diagnostic has column-basis equality

\[
R=MA,\quad R_j=\sum_i M_iA_{ij}.
\]

Let `m` and `r` denote column vectors of Matchete and RGBeta couplings. Under the **additional hypotheses** that the two programs use exactly the same physical fields and that the Matchete quartic term appears with `+m` in the Lagrangian whereas RGBeta uses `-V(r)`, equality requires

\[
M m=-R r\quad\Longrightarrow\quad r=-A^{-1}m.
\]

The script calculates and checks `A(-A^{-1})=-I` **exactly** in each of 46 full-rank blocks, retaining exact radicals. The sign `-` is explicitly a *conditional physical-convention assumption*; it is not inferred from the number `6`.

Examples (columns/rows in the order contained in each block's JSON record):

- Direct 1D block, `A=[1]`: `r=-m` **if the hypotheses hold**.
- `B_-1`/`C_-1`, corrected `110101` block:

  \[
  A=\begin{pmatrix}\tfrac12&1\\-\tfrac{\sqrt3}2&0\end{pmatrix},\qquad
  -A^{-1}=\begin{pmatrix}0&\tfrac{2\sqrt3}{3}\\-1&-\tfrac{\sqrt3}{3}\end{pmatrix}.
  \]

- `E_+0`/`E_-2` 2D triplet-only blocks:

  \[
  A=\begin{pmatrix}\tfrac12&1\\\tfrac{\sqrt{33}}6&0\end{pmatrix},\qquad
  -A^{-1}=\begin{pmatrix}0&-\tfrac{2\sqrt{33}}{11}\\-1&\tfrac{\sqrt{33}}{11}\end{pmatrix}.
  \]

All these matrices remain **candidates**; they must **not** be silently inserted into the production numerical or RGE pipeline.

## Input integrity and strict failure conditions

The tool reads the user's **fresh post-fix**:

- `output/general_t3_stage22/matchete_tensors.json`
- `output/general_t3_stage22/rgbeta_native_declarations.json`
- `output/general_t3_stage22/rgbeta_native_invariants.json`
- `output/general_t3_stage22/stage22_aligned_geometry_diagnostics.json`

It checks the aligned report's SHA-256 hashes against all three native inputs, total 16 models/50 slots/46 blocks, the four post-fix geometry classes (32 direct, 8 pseudoreal, 6 triplet-CG, 4 corrected B/C), exact invariant function agreement, registered normalisation factor, valid unique native field permutations, geometric component residual zero, invertibility, and full slot coverage.

Do **not** use a pre-fix `stage22_aligned_geometry_diagnostics.json`. The Python tests can synthesize a corrected fixture from older saved evidence solely for regression testing; **synthetic fixtures are never treated as live exports**.

## PowerShell (repository root)

Install `T3_General_Scalar_Stage23_Prefactors.zip` locally:

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage23_Prefactors.zip -DestinationPath . -Force

python -m pytest -q tests\python\test_general_t3_stage23_prefactors.py

python -m scripts.Stage23QuarticPrefactorAudit
```

The last command requires the actual corrected native exports already generated during Stage 22, including `stage22_aligned_geometry_diagnostics.json`. The script writes:

`output/general_t3_stage22/stage23_quartic_prefactor_candidates.json`

Expected message: `50/50` native prefactors, `46/46` exact conditional map blocks, **UNCERTIFIED** physical mapping. Existing Stage 22 certification remains false.

## Remaining certification gate

A physical certificate requires an independently verified common **real-scalar field-component basis**, including the `B/C` pseudoreal reidentification and triplet CG intertwiners, and an explicit live-RGBeta/Matchete scalar-action comparison fixing the **relative sign** and `PlusHc` counting. Until then:

- `source_prefactors_verified=true`: the **specified source factor/averaging identities** and exported native factors passed;
- `global_su2_component_basis_verified=false`;
- `live_real_scalar_quartic_action_verified=false`;
- `signed_physical_coupling_map_certified=false`.
