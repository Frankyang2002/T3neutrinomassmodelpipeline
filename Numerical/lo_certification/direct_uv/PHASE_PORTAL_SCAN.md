# Workstream 4: exact direct UV phase transport and native portal alternatives

## Scope and certification status

**Status: conditional numerical certification only.** This batch extends the existing 16/16 provisional exact one-loop calculation; it does not claim a physical Matchete–RGBeta–Weyl dictionary, introduce a new EFT threshold, or modify `Numerical/direct_one_loop/core.py`.

Inspected source on `main` (read only):

- `Lagrangian/interactions/T3Topology.wl` constructs `PlusHc[y1 NCM[Bar[L],CConj[F]] S1]`, `PlusHc[y2 NCM[Bar[L],F] S2^dag]`, and `PlusHc[lambdaT3 H H S1 S2^dag]`. The source also allows an opt-in independent conjugate-scalar Yukawa for hypercharge-zero S1/S2.
- `Lagrangian/interactions/ScalarPotential.wl` constructs the *native* `lambdaH{1,2}Inv{k}` portal basis as `1/2 PlusHc[...]`, not a proved equality to the benchmark `lambdaH{1,2}`.
- `output/general_observable/lo_direct_one_loop/source_audits/native_portal_identity_decomposition.json` reports 32/32 *geometric* `H†H S†S` identities representable by these source projectors, with physical coupling mapping explicitly pending.
- `Lagrangian/RunT3YukawaSpinorPhaseBridge.py` supplies a chosen explicit Weyl charge-conjugation convention for A(0), but explicitly does **not** establish the installed Matchete abstract GammaCC sign, the RGBeta epsilon component orientation, or a global Yukawa/MF phase transport.

## Equations and what is independently checked

The unmodified exact engine accepts

\[
\mathcal L_{m}=-\tfrac12\psi^{T}M_F\psi+\mathrm{h.c.},\quad
\mathcal L_Y=-\nu_\alpha\psi_i(Y_{\alpha is}\phi_s+Z_{\alpha is}\phi_s^*)+\mathrm{h.c.},
\]
\[
V_2=\phi^\dagger A\phi+\tfrac12(\phi^T B\phi+\mathrm{h.c.}).
\]

For any diagonal **unitary** field redefinition \(\psi=D_F\psi'\), \(\phi=D_S\phi'\), the transformed coefficients are

\[
M_F'=D_F^T M_FD_F,\quad A'=D_S^\dagger A D_S,\quad B'=D_S^TB D_S,
\]
\[
Y'_{\alpha is}=Y_{\alpha is}(D_F)_{ii}(D_S)_{ss},\quad
Z'_{\alpha is}=Z_{\alpha is}(D_F)_{ii}(D_S)^*_{ss}.
\]

The new `phase_portal.py` independently verifies a **fully complex action identity** (mass, scalar potential, Y and Z) and exact **signed** field-redefinition orbit(s) through the original full UV loop evaluator. General continuous complex rotations of near-cancelling matrices are not falsely presented as numerical full-loop tests because the original Takagi symmetry tolerance can be sensitive to floating-point cancellation. No existing core implementation is changed.

For each scalar \(S_i\), with native full-tensor identity-span coefficients \(c_{ik}\) and neutral-Higgs projected diagonal entries \(P_{ik,s}\), this module verifies

\[
\max_s\left|\sum_k c_{ik}P_{ik,s}-1\right|<2\times10^{-9}.
\]

For the *restricted frozen ordinary benchmark* only, its two explicitly conditional portal hypotheses are

\[
A_{ss}^{(\sigma)}=m_{S_i}^2-\sigma\,\lambda_{H_i}^{\rm RGBeta}\frac{v^2}{2}
\sum_kc_{ik}P_{ik,s},\qquad \sigma=\pm1.
\]

The \(\sigma=+1\) hypothesis **reproduces** the existing provisional direct one-loop inputs and is checked against every scalar diagonal. The \(\sigma=-1\) alternative is an **unverified** cross-package Lagrangian-sign candidate. It is propagated through the exact non-expanded loop with both scalar multiplets and all heavy fermions present. No other potential terms or source CGs are guessed.

Also calculated: (1) flipping only the T3 S1–S2 mixing sign, (2) flipping only the y2 vertex sign, and (3) their simultaneous S2 field-coordinate transformation. The latter is physically redundant; the two single-sign hypotheses are source-basis equivalent to each other when source-extra couplings are not supplied. Each conditional calculation retains vertex-pair complex matrices and Takagi masses; the alternative portal also retains the charge/heavy-fermion/scalar suppressions and cancellations via the existing `partition_exact` decomposition. All original baseline individual loop contributions remain in `models/*.json`.

**These sign alternatives are not source certified physical predictions.** In particular, neither the exact benchmark portal sign nor the independent Yukawa phases nor symmetry-allowed extra couplings is fixed by these audits. Six models (`A_-2`, `A_+0`, `D_-2`, `D_+0`, `E_-2`, `E_+0`) admit an additional conjugate-scalar Yukawa not supplied in the restricted benchmark. The direct comparison must not be forced to agree with the EFT outside expected expansion accuracy.

## Installation (PowerShell, from repository root)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$HOME\Downloads\T3_LO_direct_UV_phase_portal_batch.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
python -m pip install numpy sympy
python -m unittest discover -s tests -p 'test_lo_direct_uv_*.py' -v
```

Run against your **existing** source-derived 16 provisional inputs and native portal source audit, writing **only a new output directory**:

```powershell
python scripts/RunLODirectUVCertification.py `
  --project-root . `
  --convention-scan `
  --output-dir output/general_observable/lo_certification_direct_uv/phase_portal_01
```

The launcher exits nonzero if any of the 16 full-UV baselines or any of the 16 requested signed-portal scans is blocked. If the source audit is missing, stale or does not match the frozen benchmark SHA-256, the base UV matrix is still retained but the conditional scan is marked **blocked**, never silently substituted.

### Outputs (all new)

- `output/general_observable/lo_certification_direct_uv/phase_portal_01/all_16_certification_status.json`
- `output/general_observable/lo_certification_direct_uv/phase_portal_01/all_16_summary.csv`
- `output/general_observable/lo_certification_direct_uv/phase_portal_01/all_16_convention_scan_summary.csv`
- `output/general_observable/lo_certification_direct_uv/phase_portal_01/models/<MODEL>.json`: all existing baseline exact terms, with the new scan status and source hashes
- `output/general_observable/lo_certification_direct_uv/phase_portal_01/convention_scan/<MODEL>.json`: exact alternative matrices, per-vertex pairs, Takagi masses, conditional minus-portal suppressed-charge decomposition, source projector identities, signed-orbit residuals and explicit blockers

No previous results are overwritten. An output directory containing any model, scan or status file is refused.

## Local tests and limits

The additive suite has **32 passing tests** in an isolated source-compatible numerical harness (22 previously existing plus 10 additional). Tests exercise 16 model charge structures and the 16-model *synthetic* runner path with SHA-anchored source-shaped portal records, complex action transport, signed-loop invariance, independent Yukawa/mixing sign equivalence, and missing/mismatched identity negative controls. These are **not** 16 live Matchete source executions. Run the PowerShell command to obtain the true local-source all-16 conditional portal results.

Physical status remains **0/16** until independent native 4-spinor-to-Weyl amplitude normalization, mass phase transport, and signed portal/mixing coupling maps are furnished. The correct EFT ordering remains **UV → integrate all BSM fermions at \(\mu_F\) → EFT1 with both BSM scalars → integrate both scalars at \(\mu_S\) → SMEFT**; SM fermions are not integrated out.
