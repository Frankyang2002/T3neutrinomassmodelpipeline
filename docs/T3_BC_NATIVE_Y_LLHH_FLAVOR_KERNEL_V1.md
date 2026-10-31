# Native finite Y/LLHH: exact neutral-component flavour-spurion kernel v1

## Scope and input evidence

Consumes the **already completed** native scalar matching audits, with unchanged
source hashes:

- `output/audits/t3_bc_scalar_finite_term_extraction_v1/summary.json` and Y/LLHH
  scale-equal-mass `*.inputform.txt`, including their recorded source hashes;
- `output/audits/t3_bc_pole_log_consistency_v2/summary.json` and all **three**
  Y/LLHH quartic fragments for each model;
- the original saved native Y/LLHH Matchete expression (unless the explicit
  `--skip-original-source-validation` switch is supplied).

The v1 pole/log/constant splitting and v2 formal pole/log identity have already
been tested independently. This stage does **not** interpret the pole omission
as a renormalized MS-bar counterterm.

## Calculation actually implemented

Parse each of the 32 formal finite additive terms into exact rational monomials
with opaque Matchete fields and spinors. Require, **per term**, two Higgs fields,
two lepton fields inside one `NCM[...GammaCC...]`, two `eps[SU2L]` CG factors,
one scalar quartic, two lepton Yukawas, four mixing spurions `T3MEUij`, one
`hbar`, and one inverse `MF`. Unexpected algebra/field inventory fails closed.

Set H to component 2 and the lepton doublet to neutrino component 1. Contract
explicit two-dimensional antisymmetric tensors with epsilon(1,2)=+1.
Preserve the original `Proj[+1]` (`PR`) and `Proj[-1]` (`PL`) operators
**separately**. Treat `Bar[Coupling[...]]` as formal complex conjugation,
and the two `Flavor` slots as *ordered* indices i,j (do not silently take
`(K+K^T)/2`).

The resulting component coefficient is a sum of terms of the form

    (hbar/MF) * rational * 3^(exponent/2) * lambda^(*)
       * product(U_ij or conj(U_ij)) * y_{a,i}^(*) * y_{b,j}^(*)

with the `Proj[±1]` spinor bilinear **not** converted to a canonical Weinberg
operator. B and C have different CG and `sqrt(3)` conventions and are audited
independently.

Each model yields 16 exact spurion monomials in `PR` and 16 in `PL`, across
three quartic source labels: `lambdaT3`, `lambdaT3ExtraQ020020I1`, and
`lambdaT3ExtraQ020200I1`. The two chirality **spurion prefactors** are exact
formal conjugates under the source conventions. This is *not* a proof of
Majorana-bilinear normalization or of physical EFT matching.

## Commands (PowerShell, project root)

    Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
    python -m unittest discover -s tests -p "test_t3_bc_native_y_llhh_kernel_v1.py" -v
    python -u -X utf8 -m scripts.ProjectT3BCNativeYLLHH --all

Source and output paths are kept inside the project. No Matchete rerun, RG
rerun or change to the saved numerical neutrino masses occurs.

Outputs:

- `output/audits/t3_bc_native_y_llhh_flavor_kernel_v1/B_m1_native_y_llhh_neutral_flavor_kernel.json`
- `output/audits/t3_bc_native_y_llhh_flavor_kernel_v1/C_m1_native_y_llhh_neutral_flavor_kernel.json`
- `output/audits/t3_bc_native_y_llhh_flavor_kernel_v1/summary.csv`
- `output/audits/t3_bc_native_y_llhh_flavor_kernel_v1/summary.json`

The Python `evaluate_spurion_kernel` utility accepts *explicit* U, lambda,
flavour y1/y2, MF and symbolic-hbar inputs and returns independent PL and PR
3x3 component-coefficient matrices. It does not fetch/guess running boundary
parameters or identify hbar with 1/(16 pi^2).

## Not yet certified / next required work

1. The real Matchete one-loop renormalization convention, local counterterm
   subtraction, and matching with the EFT-side loops.
2. Canonical Weinberg-operator normalization, including the Majorana bilinear,
   chirality conjugation, symmetric flavour projection and SM component factors.
3. True running couplings/mass-mixing inputs at `muY`, including an independently
   certified quartic map from Matchete to the numerical tensor conventions.
4. Inherited Y/LLXX contribution, X-threshold loop matching and intermediate RGE;
   the earlier raw X/LLHH constant is zero **only in the source's formal
   equal-scale expression**, not a statement that physical X matching vanishes.
5. The original F-first EOMEFT equivalence remains unproven.

**No** `C5`, neutrino mass, PMNS parameter or `physical_matching_certified`
flag is modified by this stage.
