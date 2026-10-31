# Workstream C v17 — independent signed triangle master versus formal v16

Paths:

- `scripts/CompareT3BCFiniteEnvelopeToTriangle.experimental.py`
- `scripts/RunT3BCIndependentTriangle.ps1`
- `docs/WorkstreamC_IndependentTriangle_v17_README.md`

## Scope and purpose

This does **not** certify full renormalized sequential EFT matching. It performs
an **independent zero-momentum three-propagator integral** for the
`fermion F + scalar X/Y + quartic` one-loop topology and compares its strict
heavy-fermion and light-X expansion with the **source-verified** v16 formal
coefficient. It provides a specific quantitative test of finite constants,
not another term inventory.

With `A=M_F^2`, `B=M_Y^2`, `C=M_X^2`, define the signed master function

`Phi(A,B,C) = -A * integral_0^infinity [t dt / ((t+A)(t+B)(t+C))]`.

The sign is selected to compare against Matchete's signed scalar-quartic
source; this integral alone does **not** establish a complete EFT vertex
normalization. Exact analytic evaluation uses

`Phi = -A*[ B*ln(A/B)/(A-B) - C*ln(A/C)/(A-C) ]/(B-C)`

with the continuous derivative at `B=C`. The self-test checks the analytic
formula against direct positive-domain numerical integration, including equal
masses, plus the strict `1/M_F` expansion.

For `A >> B,C`:

- `Phi_XX = 1 + ln(C/A)`;
- `Phi_YY = 1 + ln(B/A)`;
- `Phi_XY = ln(B/A) - [r ln(r)/(1-r)]`, `r=C/B`.

At the **strict Y hard-boundary `r -> 0` limit**, `Phi_XY -> ln(B/A)`.
For the actual benchmark ratio `r ~= 0.76`, the omitted term is large.

Using actual Matchete PL projected three-Yukawa spurion polynomials
`P_F,P_Y,P_X,K_F,K_Y,K_X`, the script checks the nontrivial equality
`K_F=P_F`, already separately checks `K_X=0`, and requires the v15 verified
`P_F+P_Y+P_X=0`. Only after these checks does it define the **signed**
scalar pair coefficients

- `W_XX = -P_X`
- `W_XY = -K_Y`
- `W_YY = -P_Y + K_Y`

and compare the independent strict pair-master sum to **all three** ordered
Yukawa structures in the v16 formal envelope. Incorrect finite constants
cause the check to fail. The inferred pair coefficients refer only to the
previously audited scalar-quartic/LLSS topology.

The `r != 0` XY correction to the strict source is

`Delta K = K_Y * r*ln(r)/(1-r)`.

This is a **full-mass triangle diagnostic**, **not an additional matching
term to inject into v16**. Such insertion would mix EFT truncation orders,
and the complete UV/EFT subtraction, field/parameter/EOM finite shifts,
EW running, and full one-loop source equivalence are not certified.

## Commands

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WorkstreamC_IndependentTriangle_v17.zip" -DestinationPath . -Force
& .\scripts\RunT3BCIndependentTriangle.ps1
```

Main output: `output/audits/t3_bc_triangle_master_v17/summary.json`.

Optional strict/finite-X/full-mass comparison for explicit Yukawa/quartic
inputs (same schema as v16):

```powershell
python scripts\CompareT3BCFiniteEnvelopeToTriangle.experimental.py `
  --numeric-input path\to\benchmark_numeric_inputs.json `
  --output output\audits\t3_bc_triangle_master_v17\numeric.json
```

Every produced matrix is explicitly marked **NOT PHYSICAL**. Neither
production matching code nor GitHub is modified.
