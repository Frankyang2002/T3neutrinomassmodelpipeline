# Workstream A: 50-slot scalar *extra-quartic* dictionary (conditional, 9 Nov 2026)

**Input provenance:** Actual files uploaded from the user's Windows execution:

- `output/parallel_t3/wilson/matchete_native_tensors.json`: Matchete Stage22 exact CG components, 16 models / 50 extra quartic slots.
- `output/parallel_t3/wilson/rgbeta_abstract_components_conditional_v2.json`: abstract RGBeta invariants evaluated in a **chosen** symmetric-pair SU(2) frame; 16 models / 50 extra quartic slots. Its own flags say `prefactor_certified=false` and `matchete_scalar_basis_aligned=false`.

The public T3 repository `main` was read-only checked at `a944e0030b8ef24b23342d99e0d8f33cef3cffd5`. The user's uncommitted Windows source is not observable in this runtime; the hashes of the two *exported evidence files* are recorded in every result. The supplied JSON fields are treated as data, not as executable instructions. No GitHub modification is performed.

## What this batch actually calculates

The new code **does not equate the Matchete and RGBeta tensors' component indices**. It uses the *actual A_+0 Matchete CG tensor* for `Hdag Hdag S1 S2` as a **conditional** Cartesian SU(2) triplet reference. Let `M_ijA` denote that 2×2×3 tensor. From the uploaded numerical components,

\[
G_{ijA}=M_{ijA}/\sqrt{2/\sqrt3},\quad
\sum_{ij}G^*_{ijA}G_{ijB}=\delta_{AB},\quad
\sqrt{2/\sqrt3}=1.0745699318\ldots.
\]

This is **Matchete-aligned by construction**, not a native RGBeta phase or CG certificate. The source-exported symbolic RGBeta invariants (using only `del`, `eps`, `delS2`, `T3RGBetaS2TripleInvariant`) are then re-evaluated in that chosen frame. The two `delS2` factors connecting `Bar[H] Bar[H]` to `S2 S2` are conjugated on the unbarred pair; otherwise the symbolic triplet contraction is **not** the same invariant in Cartesian complex coordinates. The evaluator independently checks the component polynomial at every allowed scalar index.

For each set of identical scalar species (including bars), the code solves

\[
 R A=M, \qquad
 (\lambda_R)_j=\sum_i A_{ji}(\lambda_M)_i
\]

on the raw **complex monomial tensors**, rather than fitting a single real coefficient per coupling. The independent numerical holdout checks both the complex coefficients and their Hermitian conjugates:

\[
 V_4\big|_{\text{extra}}=2\operatorname{Re}\sum_i \lambda_i\mathcal I_i.
\]

It tests full rank, rank equality, coefficient residuals and random complex-polynomial residuals. It refuses incomplete shape, unexpected `6*Matrix[...]` native quartic metadata, corrupted source values, missing invariants, or changed conventions. No input files are overwritten.

### Results from the *uploaded* 50-slot files

| Diagnostic | Actual result |
|---|---:|
| Extra quartic slots / models | **50 / 16** |
| Maximum component-basis relative fit residual | **4.44×10⁻¹⁶** |
| Maximum independent complex-potential holdout residual | **5.37×10⁻¹⁵** |
| Extra-quartic subspace rank (each model) | Full for all nonempty models |
| Model-points with no extra quartics | A_-4, D_+2, E_-4, E_+2 |
| A(-4) portal and baseline potentials included here | **No** |

The non-diagonal 2×2 maps *in the conditional frame and unit raw-monomial normalization* are:

\[
 A_{B_{-1}}=A_{C_{-1}} =
 \begin{pmatrix}0&-2/\sqrt3\\ 1&1/\sqrt3\end{pmatrix},\qquad
 A_{E_{+0}}=A_{E_{-2}} =
 \begin{pmatrix}0&2\sqrt{3/11}\\1&-\sqrt{3/11}\end{pmatrix}.
\]

The exact source/target symbol ordering is stored in each per-model JSON group. Other groups map diagonally, with certain triplet-spurion factors \(\sqrt{2/\sqrt3}\) (A/D) and \(\sqrt2\) (E). These factors are **raw tensor** factors, not established native scalar-Lagrangian coupling conversion factors.

## Physical status — mandatory limitations

1. **Not a native RGBeta component export.** The source RGBeta package stores *symbolic* invariant tensors. Our target coordinates are deliberately inferred from the Matchete anchor, so the tiny residual is a test of *conditional algebraic compatibility*, not independent native normalization.
2. **Raw monomial vs physical quartic action.** RGBeta's saved `Coupling` function is `6*Matrix[...]` for each of the 50 records, along with `UniqueArrangements` permutation data; the original `AddQuartic` implementation uses a factor `24 × (1/√2)^4 = 6` for four complex fields. `AvgPermutations` averages arrangements. It is plausible that these restore a unit complex-field monomial potential, but a **native end-to-end contraction of the quartic real-scalar action against the Matchete Lagrangian remains unverified**. We do not multiply the fit by 6 or claim normalization is certified.
3. **Only the *extra* quartics.** Four model-points have zero extra slots. The general baseline scalar quartics (`lambdaH`, topology quartic, Higgs–scalar portals, etc.), quadratic masses, Yukawas, and their conjugation/phases have not been mapped. Crucially, A_-4 `lambdaH2Inv1` is anti-Hermitian **before multiplication by `I` inside `PlusHc`** and is not present in this 50-slot export. Do not undo or modify its verified correction in `Lagrangian/interactions/ScalarPotential.wl`.
4. **No signed GammaCC→Weyl factor.** The full one-loop dimension-five anomalous dimension is still disputed (six versus corrected seven-term). The coupled Wilson numerical ODE remains usable only as `Candidate_Qualitative`, not a physical neutrino mass prediction.
5. **EFT order unchanged:** UV → integrate out *all* heavy fermions at \(M_F\) → coupled full EFT1 Wilson evolution with *both* scalars → simultaneously integrate out *all* scalar eigenstates at \(M_S\) → SMEFT. Do not append the old fixed-order LLSS→C5 leading log to the evolved C5.

**Do not feed the v3 map automatically to Workstream B couplings or Workstream C matching as a physically certified conversion.** It can be used under an explicitly labelled conditional physics mode with matching input assumptions.

## Installation and exact PowerShell commands

Extract the one ZIP from this answer into the local repository. All four paths are new A-owned files. They do not overwrite shared Lagrangians or B/C files.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "C:\path\to\t3_workstream_A_conditional_tensor_dictionary_v3.zip"
Expand-Archive -Path $zip -DestinationPath . -Force

.\.venv\Scripts\python.exe -m unittest discover `
  -s tests -p "test_t3_parallel_a*.py" -v
```

Run on the two input files **already generated on the Windows machine**:

```powershell
.\.venv\Scripts\python.exe -m scripts.RunT3ParallelAFromNativeEvidence `
  --matchete output\parallel_t3\wilson\matchete_native_tensors.json `
  --rgbeta output\parallel_t3\wilson\rgbeta_abstract_components_conditional_v2.json `
  --output output\parallel_t3\wilson\additional_quartic_dictionary_conditional_v3.json
```

Alternatively, use the uploaded archive directly if you retained it locally:

```powershell
.\.venv\Scripts\python.exe -m scripts.RunT3ParallelAFromNativeEvidence `
  --archive output\parallel_t3\wilson\A_tensor_dictionary_inputs.zip `
  --output output\parallel_t3\wilson\additional_quartic_dictionary_conditional_v3.json
```

Each command refuses to overwrite an existing result. Use a new filename for a repeated run. Print the important summary:

```powershell
$d = Get-Content output\parallel_t3\wilson\additional_quartic_dictionary_conditional_v3.json -Raw | ConvertFrom-Json
$d.status
$d.max_component_fit_residual
$d.max_independent_potential_holdout_residual
$d.models.'B_-1'.groups | Where-Object { $_.source_indices.Count -eq 2 } |
  Select-Object species_multiset, map_rgbeta_from_matchete
```

## Next physical step

Run a native, source-consistent RGBeta **real-scalar quartic action** check (including `Coupling`, `UniqueArrangements`, complex/real scalar projectors, and H.c.) and compare against Matchete through a proven scalar basis transformation. Only after that can the conditional monomial map be promoted to a verified Matchete↔RGBeta *coupling* dictionary. Then obtain the signed GammaCC→Weyl normalization independently before promoting Wilson outputs to physical predictions.
