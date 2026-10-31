# Workstream C: Y scalar-bubble finite-structure audit (experimental v4)

## Scope and findings from the previous local run

The user-provided local output confirms that the v3 native Matchete CG export and complex-scalar tree-kernel contraction completed for **B(-1)** and **C(-1)**. For each model, 50 Y/LLHH pole coefficient structures were compared, and every ratio was `1`. The older Matchete `CheckLagrangian::Hermiticity` warnings belong to the independently component-verified original-vs-canonical representation problem; they were not newly fixed by this test. **Do not rerun the CG exporter unless you need fresh CG data.**

The v4 script tests the next _distinct_ question:

- From the original EFT1 `LLSS` and `HHSS` vertices, independently evaluate the complex-scalar Wick contractions using the v3 *native CG arrays*. Keep `YY` and `XY` kernels **separate**. Pure `XX` loops are excluded from this Y-related diagnostic.
- Check `K_YY + K_XY == P_Y`, where `P_Y` is the saved native hard Y/LLHH pole residue, stripping the shared `hbar=1/(16 pi^2)` and `1/MF` factors exactly as v3.
- Check `K_XY == F_Y(mu=M_Y)`, where `F_Y` is the saved **formal** Matchete Y finite expression (`NOT_PHYSICAL_C5`), not a renormalized Wilson coefficient. This is a new exact check, **not** assumed to pass: mismatches are saved in the output JSON and cause exit code `1`.
- Independently calculate the finite part of a generic massive zero-momentum scalar bubble as a mass-ratio **reference**. This is **not** automatically a correct replacement for the EFT matching expression.

This script does not call Wolfram, does not rewrite original files, and does not implement EFT counterterm subtraction or `GammaCC` -> symmetric Weyl Weinberg normalization.

## Scalar bubble calculation

For scalar propagators with positive real masses, define

\[
I_{ab}(\bar\mu)=\bar\mu^{2\epsilon}\int \frac{d^{4-2\epsilon}k}{(2\pi)^{4-2\epsilon}}\frac{i}{(k^2-m_a^2+i0)(k^2-m_b^2+i0)}.
\]

Up to the common Feynman-sign convention, the **dimensionless** Laurent kernel is

\[
16\pi^2 I_{ab}\sim \frac{1}{\bar\epsilon}
 + \int_0^1 dx\log\frac{\bar\mu^2}{x m_a^2+(1-x)m_b^2}.
\]

At \(\bar\mu=M_Y\), for \(m_a=m_b=M_Y\), the finite kernel is `F_YY=0`. For \(m_a=M_X,m_b=M_Y\), set \(r=M_X^2/M_Y^2\). The full-mass bubble has

\[
F_{XY}(r;M_Y)=1+\frac{r\ln r}{1-r},\qquad 0<r<1.
\]

The strict `mX -> 0` limit gives `F_XY=1`, which motivates comparing the saved constant finite terms with `K_XY`. For the B/C benchmarks \(r\simeq0.758\)-0.759, the full massive bubble reference instead has `F_XY ~ 0.13`. **That difference is not a derived correction to physical `C5`.** A finite-mass sequential-EFT matching prescription must consistently include hard/soft regions, EOM/higher-dimensional operators and the EFT2 subtraction. A near-degenerate scalar spectrum makes truncating terms in `MX^2/MY^2` numerically consequential.

The output's `scalar_bubble_msbar_finite_at_muY` values must not be multiplied into `K_XY` and relabeled a renormalized Wilson coefficient.

## Installation

Extract `T3_WorkstreamC_YFiniteScalarBubble_v4.zip` into the repository root. It creates:

- `scripts/CheckT3BCYFiniteScalarBubble.experimental.py` (complete new file).
- `docs/WorkstreamC_YFiniteScalarBubble_README.md` (this file).

It **requires** the previous v2/v3 scripts and the already created local Y/CG audit files. It does not replace them.

## Exact PowerShell commands

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "$HOME\Downloads\T3_WorkstreamC_YFiniteScalarBubble_v4.zip" `
  -DestinationPath . -Force

python scripts\CheckT3BCYFiniteScalarBubble.experimental.py --self-test

python scripts\CheckT3BCYFiniteScalarBubble.experimental.py `
  --output output\audits\t3_bc_y_finite_scalar_bubble_v1\summary.json
```

Required existing project-relative input files:

- `scripts/CompareT3BCYNativeSU2Coefficients.experimental.py` (v3)
- `scripts/AuditT3BCYMixedScalarPoleSupport.experimental.py` (v2)
- `output/audits/t3_bc_y_native_cg_components_v1/{B_m1,C_m1}_native_CG.json`
- `output/audits/t3_bc_y_native_cg_components_v1/comparison.json`
- `output/audits/t3_bc_y_renormalization_ledger_v1/{B_m1,C_m1}_y_llhh_raw_pole_residue.inputform.txt`
- `output/audits/t3_bc_y_renormalization_ledger_v1/{B_m1,C_m1}_y_llhh_formal_finite_at_muY.NOT_PHYSICAL_C5.inputform.txt`
- `output/general_t3_mass_eigen_matchete/eft1_fermion_first_tree/{B_m1,C_m1}_final_eft1_fermion_first_tree_llss_tree.inputform.txt`
- `output/general_t3_mass_eigen_matchete/eft1_fermion_first_tree/{B_m1,C_m1}_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt`
- `output/general_t3_mass_eigen_matchete/{B_m1,C_m1}_final_parameters.json`

If the previous Y ledger is absent locally, it must be built with your already installed ledger script and original `xy_fresh_native_evidence_v3.zip`; no native matching need be rerun.

## Status boundary

- Source-level pole-to-tree contraction: **user v3 test passed**, both models, ratio `1`.
- Tree-to-formal-finite XY coefficient equality: **not verified until v4 is run locally**.
- Full-mass bubble finite reference: **analytic**, checked by independent numerical Feynman-parameter quadrature in v4 self-test.
- EFT1/EFT2 independent operator counterterms, UV/IR separation, finite threshold scheme, `GammaCC` normalization, running, X-threshold inheritance and physical `C5`: **not computed/certified**.

Next after v4: independently derive and verify the EFT1/EFT2 UV operator counterterm **difference** (including any soft-EFT overlap) before converting a finite scalar matching expression into physical `C5`.
