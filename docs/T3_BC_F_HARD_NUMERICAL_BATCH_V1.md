# F-hard finite-source numerical batch (B(-1), C(-1))

## Scope
UV -> EFT1 -> SM, one intermediate EFT. No production code modifications.
For each of three running heavy Majorana Takagi eigenstates, evaluate all 68
pole, 68 mass-log and 68 constant native Matchete F-hard LLHH terms, using
verified native neutral component reduction, complex conjugation and source
mixing spurions. Inputs are from the SAME five-mass UV state at muF = 1e5 GeV.

Let \(L_r=(\bar y_{1r},\bar y_{2r})\), and let \(K_s(U,\lambda)\) be the
source-derived ordered 2x2 coefficients in sector s. The independent PL
ordered matrix is \(A_{s,r}=L_r K_s L_r^T/(16\pi^2 M_r)\) and
\(C_{s,r}=A_{s,r}+A_{s,r}^T\). Do **not** add PR again. The *formal*
source quantity is
\[
C_F^{\mathrm{formal}}(\mu_F)=\sum_r\left[C_{\mathrm{constant},r}
+\ln(\mu_F^2/M_r^2)C_{\mathrm{log},r}\right].
\]
The exact source audit establishes C_pole=C_log, and the numerical derivative
checks \(dC_F^{\rm formal}/d\ln\mu_F=2\sum_r C_{\rm pole,r}\)
**at fixed running endpoint inputs** (not an RG-scale independence test).

This source expression has not had the required UV/IR pole decomposition,
independent EFT1 counterterms, field/EOM finite shifts, or operator inheritance.
It is not an independently renormalised Wilson matching condition. It is
intentionally NOT injected into the previously calculated F threshold nor
added to the scalar-bubble contribution (double-counting risk).

## Run
```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\t3_bc_f_hard_numerical_batch_v1.zip" -DestinationPath . -Force
& .\scripts\RunT3BCFHardNumericalBatch.ps1
```

## Sources and outputs
Requires original source `output/general_t3_mass_eigen_matchete/eft1_fermion_first_oneloop/` and
Batch 1 `output/studies/t3_bc_fermion_threshold_batch_v1/{summary.json,f_hard_formal_source.json}`
as well as five-mass `output/studies/t3_full_stack_overnight_v1/bc_full_mass`.
Checks exact raw source SHA, exact source polynomial symbolic equality,
source chirality content, native PL/PR hermitian conjugacy, and mass-log derivative.
Writes per-model full matrices and scan to
`output/studies/t3_bc_f_hard_numerical_v1/`.

**No physical C5 certification or UV-to-EFT1 renormalized finite subtraction.**
