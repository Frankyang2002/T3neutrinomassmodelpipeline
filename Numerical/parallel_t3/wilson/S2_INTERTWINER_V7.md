# Workstream A v7 — native RGBeta S2 frame/intertwiner audit (schema fix)

## Verified source boundary

Reviewed GitHub `Frankyang2002/T3neutrinomassmodelpipeline`, `main`
`2bc2e4b516776212645370105cc63495a3c47348` (2026-11-09).
The current `T3RGBetaModel.wl` registers `SU2L[S2]` and an adjoint portal
` tGen[SU2L[fund],adj,#1,#2] tGen[SU2L[S2],adj,#3,#4]`.
`rgbeta_fundamental_fierz_native_v6.json` records the **symbolic** `tGen`
for S2; it does not contain three numerical S2 generators in specified axes.
The archived B2 native RGBeta `orientation_witness.json` says
`NATIVE_SYMBOLIC_GENERATOR_RULES_CAPTURED_NUMERIC_AXES_PENDING`.
It is a *historical diagnostic*; its old source HEAD/hashes are not presented
as a fresh v7 native execution.

**Important source-provenance trap:**
`output/general_observable/lo_direct_one_loop/source_audits/generator_candidate_probe.json`
sets `generator_basis_is_live_export=false` and identifies the input as repository
**test fixtures**. Its downstream `matchete_generator_analysis.json` is not
independent live generator data. The v7 code must never promote it.

## v7 schema correction

The first v7 archive incorrectly read the **v5** triplet key
`candidate_rgbeta_from_native_matchete_real_matrix` from a **v6** sector.
The correct v6 triplet key is `candidate_rgbeta_from_matchete_real_matrix`.
Doublet and singlet sectors retain `rgbeta_from_matchete_real_matrix`.
The fix covers both direct sector audits and saved-v6 source revalidation,
adds explicit sector-count/identity checks and regression tests.
No physics convention or previous evidence is changed.

## What the numerical audit actually proves

Use the *declared* orthonormal symmetric-square frame

```
C[1,1,1] = 1
C[2,1,2] = C[2,2,1] = 1/sqrt(2)
C[3,2,2] = 1
```

and `t_f^a = sigma^a/2`. Compute `G_C^a = C (t_f^a tensor 1 + 1 tensor t_f^a) C†`
(in indices, with the canonical conjugations). Test isometry, symmetric
projector, equivariance, Hermiticity, trace index `T(S2)=2`, Casimir `C2=2`,
and `[G_C^1,G_C^2]=i G_C^3`.

For each of the 14 v6 triplet portal sectors, reconstruct `K_M` from **native
Matchete CGs** and the source-audited v6 conditional map. Extract
`G_M^a = 2 tr_H(t_H^a K_M)`. Check the triplet algebra and solve
`G_M^a U = U G_C^a` by the *generator intertwiner nullspace* (no potential
fit). Check `U† U = I` and transport `sum_a t_H^a tensor G_C^a` to the
Matchete CG tensor. Revalidate the v5 and v6 source hashes and component
claims before accepting any comparison.

**This `U` is NOT an independent RGBeta-to-Matchete bridge.** `G_M` is
reconstructed from Matchete's portal CG itself, and `G_C` is a user-declared
mathematical `Sym²(fund)` frame, not native measured RGBeta components.
The transport check is a valuable representation consistency check, not an
independent portal-coupling normalization certificate.

Next, explicitly rotate the chosen S2 frame by a nontrivial unitary `W`:

```
G'_C = W G_C W†
C' = W C
U' = U W†
```

Both frames reproduce all abstract S2 invariant identities and the same
Matchete tensor. Hence the symbolic RGBeta structure alone **cannot select
one unique numerical triplet frame**. It also does not independently fix
an RGBeta-to-Matchete component intertwiner, even if a chosen comparison `U`
exists and all numerical residuals are tiny.

## Files (add only)

- `Numerical/parallel_t3/wilson/s2_intertwiner_v7.py`
- `Numerical/parallel_t3/wilson/S2_INTERTWINER_V7.md`
- `scripts/RunT3ParallelAS2IntertwinerV7.py`
- `tests/test_t3_parallel_a_s2_intertwiner_v7.py`

New output (refuses overwrite):

- `output/parallel_t3/wilson/s2_intertwiner_conditional_v7.json`

## PowerShell from repository root

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = "$env:USERPROFILE\Downloads\t3_workstream_A_s2_intertwiner_v7_schemafix.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
.\.venv\Scripts\python.exe -m unittest discover -s tests -p 'test_t3_parallel_a_s2_intertwiner_v7.py' -v
if ($LASTEXITCODE -ne 0) { throw 'Workstream A v7 tests failed' }
.\.venv\Scripts\python.exe -m scripts.RunT3ParallelAS2IntertwinerV7 `
    --project-root . `
    --native-rgbeta-witness 'output\general_observable\physical_normalisation\b2_rgbeta_generator_native_02\orientation_witness.json' `
    --generator-candidate 'output\general_observable\lo_direct_one_loop\source_audits\generator_candidate_probe.json'
if ($LASTEXITCODE -ne 0) { throw 'Workstream A v7 audit blocked; inspect error' }
```

The native RGBeta and fixture inputs are *optional* historical diagnostics;
remove their command-line flags if absent. In either case the conclusion is
**conditional**. If v5/v6 evidence is stale against the current local sources,
the script **blocks** rather than silently reusing it. The v7 script does not
rerun Wolfram, scalar matching, or physical Wilson evolution.

## Next genuine certification gate

Obtain a new, independently provenance-pinned native RGBeta component basis
for `delS2` **and** `tGen[SU2L[S2]]`, with a documented numerical index
interpretation and consistent transformation of `S_i`, `Bar[S_i]` and Higgs
axes. Separately obtain an independent native Matchete generator basis (not
the test-fixture generator probe). Only then solve the *cross-package* S2
intertwiner and test held-out native invariant tensors and full normalized
scalar actions. Do not pick the RGBeta frame to fit Matchete portals.

Outstanding blockers remain: B(-1)/C(-1) mixed scalar mass diagonalization,
full scalar potential/physical coupling normalization, signed GammaCC -> Weyl,
and Wilson ADM/physical boundary certification. No physical Wilson rerun.

EFT ordering: **UV → integrate out ALL fermions → integrate out ALL scalars → SMEFT**.
