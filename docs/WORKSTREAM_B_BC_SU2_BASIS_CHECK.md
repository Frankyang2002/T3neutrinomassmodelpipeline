# Workstream B — B(-1)/C(-1) SU(2) component-basis audit

## Scope and source evidence

Read-only source inspection of the current project `main` as of 2026-10-31,
including `Lagrangian/model/T3Fields.wl`,
`Lagrangian/interactions/GeneralScalarOperators.wl`,
`Numerical/experimental/GeneralT3MixedScalarMassBasis.py`,
`RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl`,
`scripts/ProbeT3BCMixedBilinear{Matchete,RGBeta}.wl`, and
`scripts/ExportGeneralT3SU2NativeCovariance.wl`.

The installed Matchete's **live generator matrices** are acquired via its documented
`Generators[SU[2],{1}]` accessor, after loading the SM and registering the B/C fields.
The actual Matchete mixed-bilinear CG is **reused** from the earlier live
`GetCGTensor` report. RGBeta's live scalar-mass records and massive `Lam` evidence
are **reused** from the successful v5 output. The already-passed v6 conditional
cross-tool map is checked again, including live-report SHA256 provenance and
current-source digests.

No production RGE, matching, or model definition is modified.

## What is tested

Write `E = [[0,1],[-1,0]]`. Let `T_M^a` be the three native Matchete SU(2)
fundamental generators, `C_M = k E` the native registered mixed scalar bilinear
CG tensor, and `P^a = sigma^a/2` the canonical Pauli realization of RGBeta's
*abstract* SU(2) fundamental generators.

We construct the real `3 x 3` adjoint-frame rotation

    O[a,b] = 2 trace(T_M[a] P[b]),    T_M[a] = sum_b O[a,b] P[b].

The independent verifier checks Hermiticity, tracelessness,
`Tr(T_M[a] T_M[b]) = delta[a,b]/2`, orthogonality of `O`, `det(O)=±1`, the
corresponding commutation algebra, and

    T_M[a]^T C_M + C_M T_M[a] = 0,
    U_M^T C_M U_M = C_M

for 32 random finite gauge rotations. Gauge invariance of the complex scalar
bilinear and its conditional Matchete/RGBeta action is checked on random fields.
Negative tests reject wrong generator factors and non-epsilon CG tensors.

This is independent evidence of **SU(2) covariance and a consistent choice of
canonical adjoint frame**. It is not a native numeric RGBeta generator extraction:
RGBeta's stored `tGen` representation is symbolic, and the Pauli realization is
chosen explicitly for component evaluation. Nor does it identify an overall
phase for each scalar field across packages.

## Fundamental scalar phase ambiguity

If scalar components in Matchete and RGBeta are related by independent phases,

    S1_M = exp(i alpha1) S1_R,
    S2_M = exp(i alpha2) S2_R,

then, for the same holomorphic `Bar[S1] Bar[S2]` field word,

    mu_R = -k exp[-i(alpha1+alpha2)] mu_M.

The already-established `mu_R = -mu_M` corresponds to the convention
`k=1` and `alpha1+alpha2=0` (mod 2*pi).

An overall complex scalar-field phase leaves the gauge generators and kinetic
term unchanged; therefore **no experiment comparing only SU(2) generators and
an epsilon invariant can determine that relative field phase**. A universal
physical parameter-map certificate needs a complete choice of field phases,
including the Yukawa and quartic terms. Rephasing all couplings consistently
does not change physical observables. The phase-rotated action is checked
explicitly; holding the couplings fixed under an independent field rephasing
is a failing negative control.

## Install and run (PowerShell)

Repository root:

    C:\thesistemp\T3neutrinomassmodelpipeline

Extract the ZIP into that directory (substitute the download's actual path):

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive `
  -LiteralPath 'C:\path\to\Workstream_B_BC_SU2_Basis_Audit.zip' `
  -DestinationPath . -Force
python -X utf8 -m pytest -q tests\test_t3_bc_su2_basis.py
```

Generate the **new** native Matchete evidence. Choose a fresh JSON file name
if the output already exists; the tool never overwrites old evidence:

```powershell
wolframscript -file scripts\ProbeT3BCSU2GeneratorsMatchete.wl `
  output\audits\workstream_b_bc_su2_matchete_v1.json
Write-Host "Wolfram exit code: $LASTEXITCODE"
```

Use the already-successful `matchete_v1.json`, `rgbeta_v5.json` and
`mixed_bilinear_map_v6.json` from the previous Workstream B step:

```powershell
python -X utf8 -m scripts.VerifyT3BCSU2Basis `
  --generators output\audits\workstream_b_bc_su2_matchete_v1.json `
  --matchete output\audits\workstream_b_bc_mixed_matchete_v1.json `
  --rgbeta output\audits\workstream_b_bc_mixed_rgbeta_v5.json `
  --bilinear output\audits\workstream_b_bc_mixed_bilinear_map_v6.json `
  --output output\audits\workstream_b_bc_su2_basis_v1.json
```

For a compact summary:

```powershell
python -c "import json; x=json.load(open('output/audits/workstream_b_bc_su2_basis_v1.json')); print(x['status']); [print(v['model'], v['generators']['max_su2_commutator_residual'], v['invariant']['max_finite_su2_covariance_residual']) for v in x['models']]"
```

## Status / limits

- Python algebra and synthetic-input tests were executed locally.
- The new Wolfram probe **could not be run locally** (no installed WolframScript).
- Do not interpret a passing Pauli/generator compatibility check as an
  empirical determination of independent field phases or native numeric
  RGBeta generators.
- Complete renormalizable scalar-mass RG closure also needs the Higgs
  quadratic-mass beta function, a separate next step.
- Workstream A final operator-basis certificate remains independent.
