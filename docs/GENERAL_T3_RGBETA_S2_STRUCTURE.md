# General T3: Native RGBeta SU(2) `S2` and triplet conjugation structure

**Scope:** Read-only, fail-closed live RGBeta group-structure probe, plus exact Matchete triplet-metric audit. All output remains **UNCERTIFIED** for physical inter-package component-basis identification and quartic coupling normalisation. No production edits and no GitHub edits.

## Why a numerical RGBeta triplet matrix is not directly exportable

Inspected RGBeta v1.2.0 `Package/GroupsAndIndices.m` defines `tGen[rep,A,a,b]` as an **abstract representation generator** and `delS2[group,a,i,j]` as an invariant between a symmetric-square index and two fundamental indices. `RefineGroupStructures` expresses the `S2` generator in terms of fundamental **symbolic** generators. The implementation does *not* select numeric Pauli generator matrices or three explicit triplet component coordinates. `Package/FieldsAndCouplings.m` associates complex fields with a scalar/conjugate-scalar index pair, not with explicit Cartesian real/imaginary components.

A numeric triplet generator constructed from Pauli matrices and a conventional symmetric-product basis is therefore a **chosen reference representation**, not an independent native RGBeta generator-component measurement. The existence of a unitary equivalence alone does not select a unique field-phase or barred-field convention. Do not label that equivalence a physical coupling certificate.

## Independent native checks made by the Wolfram script

For `A_+0`, `A_-2`, `D_+0`, `D_-2`, `E_+0`, `E_-2`, after rebuilding each registered model, the exporter requires:

- Triplets use the actual registered `SU2L[S2]` gauge representation.
- Native dimension `Dim[S2] = 3`, quadratic Casimir `Casimir2[S2] = 2`, and trace normalisation `TraceNormalization[S2] = 2`.
- The group-level conjugate rule `tGen[Bar[rep],A,i,j] = -tGen[rep,A,j,i]` **as evaluated in the installed RGBeta runtime**.
- Native symmetric-doublet orthogonality `delS2[SU2L,a,i,j] delS2[SU2L,b,i,j] = del[SU2L[S2],a,b]`.
- The native complex-triplet/conjugate-triplet scalar contraction metric has both corresponding entries equal to 2, consistent with Stage 25.
- The native generator remains **symbolic** rather than a numerical 3×3 component matrix.

Any mismatch aborts and no JSON evidence is produced. No synthetic result may substitute for the user's actual Wolfram export.

## Exact Matchete triplet invariant form

From Stage 27's independently extracted normalised CG matrices \(C_a\) and \(\epsilon_{12}=1\), form

\[
\Omega_{ab}=\sum_{i,j,k,l} C_a^{ij}C_b^{kl}\,\epsilon_{ik}\epsilon_{jl}.
\]

For the project native Matchete CG data, **\(\Omega=-I_3\)**. The Stage 28 Python audit recomputes and checks the exact matrix and its invariance:

\[
(T^a_{\mathrm{Mat}})^T\Omega+\Omega T^a_{\mathrm{Mat}}=0,
\qquad \sum_a(T^a_{\mathrm{Mat}})^2=2I_3.
\]

This is a consistency check on the Matchete CG and a *source-compatible* abstract `S2` group structure. It does not determine the Matchete/RGBeta triplet **component coordinates**, nor does it prove the field conjugation map. The previously reported \(UU^T=\mathrm{diag}(1,1,-1)\) is not by itself the physical `Bar` map.

## Local commands

From the repository root, after extracting the ZIP:

```powershell
python -m pytest -q tests\python\test_general_t3_triplet_native_structure.py

wolframscript -file .\scripts\ExportGeneralT3NativeTripletStructure.wl `
    output/general_t3_stage22/stage28_live_rgbeta_triplet_structure.json
if ($LASTEXITCODE -ne 0) { throw "Stage28 native probe failed" }

python -m scripts.RGBetaTripletNativeStructureAudit
python -m pytest -q tests\python\test_general_t3_triplet_native_structure.py
```

Native input `output/general_t3_stage22/stage28_live_rgbeta_triplet_structure.json`; Python diagnostic output `output/general_t3_stage22/stage28_triplet_native_structure_audit.json`. The Python audit re-runs Stage 27 and Stage 26 checks before trusting the saved witness.

## Practical next step

There are two distinct goals: (a) show the RGBeta native `S2` algebra matches spin 1; this stage can verify it, and (b) fix a concrete component coordinate and physical field conjugation convention across Matchete and RGBeta; the abstract RGBeta group interface **cannot uniquely certify (b)**. To reach rough numerical RGE benchmarks efficiently, the next implementation should **explicitly declare** one shared component convention (rather than claim it was extracted from RGBeta), calculate corresponding couplings with tracked convention choice, and perform a basic basis-independent benchmark cross-check. A physically certified map would still require an independent cross-package action/field check.
