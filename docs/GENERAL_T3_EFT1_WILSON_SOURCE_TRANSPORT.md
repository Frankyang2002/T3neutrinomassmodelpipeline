# General T3: numerical extra-quartic EFT1 LLHH Wilson-source transport

## Verified preceding steps

- Stage 16: grouped two-threshold **renormalisable** UV/EFT1 one-loop RGBeta running and pre-existing baseline Weinberg evolution (both B/C(-1) run successfully on Windows).
- Stage 17: two additional scalar quartic invariants contracted into the Matchete PL EFT1 Wilson seed, producing four nonzero Higgs-real components of the **extra** scalar-pair Wilson anomalous dimension in both models (10/10 local tests and both model exports passed).
- Stage 15: potential/Lagrangian sign relation `lambda_Matchete = -lambda_RGBeta` for these couplings in the common scalar field basis. Stage18 uses the *RGBeta potential* coefficients, so it does NOT insert an additional minus sign.

## Stage18 calculation

Reuse `output/T3_{B,C}_alpha_m1/data/{uv_rgbeta_rge.json,eft1_rgbeta_rge.json}` and the genuine `output/general_t3_stage17/{B_m1,C_m1}_extra_eft1_scalar_pair_beta.json`. Run the Stage16 grouped benchmark again in memory to obtain the time-dependent EFT1 state at 17 logarithmically spaced scales. Use the full existing coupled RGEs for **all** renormalisable EFT1 couplings (including six special scalar quartics); only the Wilson insertion is new.

Stage17's scalar-pair beta components are homogeneous of degree 2 in the one-generation PL Yukawa and inverse-linear in MF. The code validates that they factor into

    beta_ijab^(1)(mu) = A_ab(mu) * conjugate(y1)^2 / MF
                       + B_ab(mu) * conjugate(y2)^2 / MF .

The only allowed free symbols in `A` and `B` are the two active complex RGBeta potential couplings and their conjugates, each to degree at most one. In the positive diagonal heavy-mass basis at the F boundary, promote the factors to symmetric 3x3 matrices:

    Y_i[p,q] = sum_r conjugate(y_i[p,r]) conjugate(y_i[q,r]) / M_r .

The Yukawas, positive heavy masses, and tree Wilson seeds are evaluated **at mu_F** and kept fixed below F, where those heavy fields are absent. Evolve the source only:

    Delta C_ijab|extra-source = (1 / (16*pi^2)) * integral_{
      ln(mu_F)}^{ln(mu_S)} beta_ijab^(1)(mu) dln(mu).

Integration uses composite Simpson quadrature over the saved EFT1 trajectory, plus an independent composite trapezoid estimate (reported as a quadrature difference). Both are downwards in logarithmic scale. Output units: `GeV^-1` for each flavour-matrix component.

### Scope (important)

This **is not the full coupled EFT1 Wilson system**: it omits homogeneous Wilson mixing (except the direct new scalar-pair source), additional operator-generated source feedback, one-loop corrections to the F threshold seed, and finite scalar-threshold matching. It exports **real** Higgs-index Wilson coefficients (`i,j` fixed to lepton SU2 indices `(1,2)`, `a,b` in Higgs real indices 1..4), not a canonically normalised Weinberg coefficient. Therefore it makes *no* new neutrino mass prediction, and does not add a contribution to the pre-existing physical baseline. This is a full **three-lepton-flavour source lift** of a tree-seed single-insertion, not a validated three-generation coupled Wilson ODE.

For radiatively generated extra quartics initialized at zero at the UV scale, this piece has nominal two-loop suppression. The computation is a selected `O(hbar^2)` effect, **not a complete two-loop prediction**.

We retain the Stage16 grouped mass hierarchy and its finite-order approximation, not individual particle thresholds. The `mu_F` tree boundary is continuous for renormalisable parameters; no new finite matching terms are assumed. The original Stage11-15 toy spectrum violates the grouped F-first ordering and is not used.

## Commands (PowerShell from project root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch18.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_e_f_t1_wilson_source_transport.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3EFT1WilsonSourceTransport.ps1
```

Optional numerical consistency scan, uses more sampled points:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3EFT1WilsonSourceTransport.ps1 -Samples 33
```

## Output files

- `output/general_t3_stage18/B_m1_extra_eft1_wilson_source_integral.json`
- `output/general_t3_stage18/C_m1_extra_eft1_wilson_source_integral.json`
- `output/general_t3_stage18/stage18_summary.json`

The script fails closed and records each model error in `stage18_summary.json`. It does not modify production physics code, the baseline outputs, or GitHub.

## Next physics task

Match these real Higgs Wilson tensors onto the correctly normalised Weinberg operator using the project's established `WeinbergTensorAdapter` and the scalar threshold, then compare with direct Matchete matching. Only after consistency checks should the extra piece be included in the final SM+Weinberg transport.
