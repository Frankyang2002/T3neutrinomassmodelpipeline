# Workstream C: independently calculated Y-quartic mixing (experimental)

**Scope:** A calculable part of the missing B(-1)/C(-1) scalar-threshold renormalization. This is deliberately **not** the full Y counterterm, a subtraction prescription, a B/C numerical result, or a physical Weinberg coefficient.

## Verified existing implementation

Read-only inspection of `Lagrangian/RunMatching.wl` (GitHub blob `d1396d4`) and `Lagrangian/RunT3BCYOnlyCanonicalReconstruction.experimental.wl` (`57afcd8`) finds that the Y script computes the one-loop-minus-tree `MatchedEFT` difference and exports the `GammaCC`/LLHH sector. It does **not** compute the independent EFT1/EFT2 counterterm difference or establish the UV/IR origin of all displayed poles. Existing `scripts/BuildT3BCYRenormalizationLedger.py` partitions and projects the native poles but also does not calculate counterterms.

The source Y outputs and archive `xy_fresh_native_evidence_v3.zip` were *not available in this session*. The audit inventory of 50 pole, 50 logarithm, and 32 constant terms for both B/C is taken from the supplied 31 October handover, not rerun here. No B/C coefficient-by-coefficient comparison was possible.

Matchete uses a hard-region functional expansion for one-loop matching; the hard/soft subtraction must not be applied twice. Artificial infrared singularities of the hard expansion can be associated with EFT ultraviolet counterterms, so identifying a pole in the raw Matchete expression does **not** establish its origin or whether it should be removed as a UV counterterm of this specific EFT.

## Independent physics result

Aebischer, Bresciani and Selimović, *Anomalous Dimension of a General Effective Gauge Theory II: Fermionic Sector*, arXiv:2512.16890v1, Eqs. (3.3) and (4.85), give the scalar-quartic contribution to the one-loop beta function of the fermion-bilinear/two-scalar dimension-five tensor.

For **real** scalar components `phi_a`, with fully symmetric quartic tensor and the paper's Weyl operator convention:

```
V_4 = lambda_abcd phi_a phi_b phi_c phi_d / 4!
L_5 = C_ijab psi_i psi_j phi_a phi_b + h.c.
dot(C)_ijab := (16*pi^2)*d(C_ijab)/d(ln(mu))

 dot(C)_ijab | quartic, internal Y
       = sum_(c,d are REAL Y components) lambda_abcd * C_ijcd
```

This follows by isolating the `+ lambda_abcd` term in Eq. (4.85). The sum is over **ordered** `(c,d)`; off-diagonal symmetric Y pairs contribute twice. The script evaluates this contraction independently using exact symbolic algebra. Its result is the quartic-induced EFT1 operator-mixing contribution, **not** the whole EFT1 beta function, EFT1−EFT2 counterterm difference, or matching coefficient.

All EFT1 scalar multiplets must first be decomposed into independent real scalar components. `lambda_abcd` must use the precise `V4/4!` convention. Native Matchete `GammaCC`/`Proj` bilinears and its Wilson coefficients must be converted into the paper's **Weyl and flavour conventions** with independently proven signs, field normalizations, combinatorial factors, and Clebsch–Gordan contractions. The conversion is **not implemented here**. Couplings may remain formal symbols, but symbolic entries must already be in the specified physical basis.

## Included files, project-relative paths

- `scripts/ComputeT3BCYQuarticMixing.experimental.py`: the full independent-contraction calculator.
- `docs/T3BCYQuarticMixing_TOY_NOT_BC.json`: a deliberately synthetic exact-arithmetic test case. **Do not use this file as B/C benchmark data.**
- `docs/WorkstreamC_YQuarticMixing_README.md`: the present instructions.

The existing production pipeline, Y matching script and `NOT_A_COUNTERTERM` evidence files are untouched.

## Reproducible PowerShell commands

First extract the supplied ZIP into the **local** project root (change the download path to where you saved it):

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WorkstreamC_YQuarticMixing_v1.zip" -DestinationPath . -Force
python -m pip show sympy
python scripts\ComputeT3BCYQuarticMixing.experimental.py --self-test
python scripts\ComputeT3BCYQuarticMixing.experimental.py `
  --input docs\T3BCYQuarticMixing_TOY_NOT_BC.json `
  --output output\audits\t3_bc_y_quartic_mixing_v1\TOY_result.json
```

The toy run should report

```
(16*pi^2)*beta_y = c1*lam1/2 + c2*lam2/2
physical_C5_Y = null
```

When and **only when** you have independently converted the original B/C `EOMEFT` LLSS coefficients and scalar quartic CG tensors into the stated real-component normalization, supply project-local input files such as

```
output/audits/t3_bc_y_quartic_mixing_v1/B_m1_real_scalar_tensors.json
output/audits/t3_bc_y_quartic_mixing_v1/C_m1_real_scalar_tensors.json
```

Then calculate, for example:

```powershell
python scripts\ComputeT3BCYQuarticMixing.experimental.py `
  --input output\audits\t3_bc_y_quartic_mixing_v1\B_m1_real_scalar_tensors.json `
  --output output\audits\t3_bc_y_quartic_mixing_v1\B_m1_quartic_beta.json
```

These benchmark input files **are not included and have not yet been produced**. Running the benchmark command before generating them will fail, appropriately.

## Input schema

`format` must be `"t3bc-y-scalar-quartic-beta-v1"`. `scalar_labels` is the ordered list of **real scalar** degrees of freedom, each with a zero-based integer index. `fermion_labels` gives the fermion flavour labels. `higgs_scalar_indices` and `heavy_y_scalar_indices` are disjoint integer arrays. The `quartic_tensor` is a list of fully symmetric entries with four `scalars` indices and a `terms` list. The `wilson_tensor` contains `flavours` and two `scalars` indices **both in Y**. Each `terms` entry has an exact `rational` string and a `factors` dictionary of formal symbol names and integer powers. `ImaginaryUnit` and `Sqrt3` are special exact factors. Unlisted tensors are zero; duplicate entries in the symmetry-reduced basis are rejected.

For illustration:

```json
{"scalars":[0,0,1,2],"terms":[{"rational":"1/3","factors":{"lam2":1}}]}
```

An optional `native_LLHH_pole_without_hbar` array can be added **only if** an independent map has converted the native LLHH coefficient into *exactly the same real-scalar, flavour, chirality and Weyl operator conventions* and has removed its explicit `hbar=1/(16*pi^2)` factor. The diagnostic compares `2P` to both signs of `dot(C)|_quartic,Y`; it does not presume either equality, and equality cannot certify counterterms, UV/IR pole assignment or physical matching.

## What is still required

1. Obtain native B/C archives and verify exact SHA-256 provenance of their source parameter snapshots, tree `EOMEFT`, CG arrays and Y matching stage diagnostics.
2. Establish the **native Matchete `GammaCC` normalization** and convert all required EFT1 `LLYY`/`HHYY` operators to the chosen real scalar/Weyl normalization. Compute the two model-specific contractions.
3. Calculate the **other operator-mixing and field-renormalization terms that survive** at this perturbative order; establish the EFT1/EFT2 scheme and their actual counterterm difference, including any EOM/finite field redefinitions required by the chosen matching convention.
4. Separate genuine UV poles from hard-region artificial IR poles using an independently specified scheme, and compare the complete counterterm-plus-hard amplitude with the native residue in every flavour/chirality/SU(2) sector.
5. Only if these pass, evaluate physical finite Y matching and certify `C5` after the native-to-Weinberg conversion and sequential inheritance.

At present, `independent_y_MSbar_counterterm_certified`, `raw_Matchete_pole_UV_IR_origin_assigned`, `net_pole_cancellation_certified` remain `false`; `physical_C5_Y` is `null`.
