# B(-1)/C(-1): running mixed quadratic, original-basis triangle and native source audit (v2)

## Purpose

This batch directly addresses the unresolved **independent complex scalar mixed quadratic** in the existing B/C running-input one-scale calculation. It checks the mass-dependent one-loop scalar kernel in the original interaction basis against an independent spectral (mass-eigenbasis) calculation, and re-rotates the actual **seven Matchete-native extra quartics** at the *running* scalar mixing matrix. The source-matched one-loop weights from Workstream C v22 are also used to quantify the numerical impact of the leading-\(1/M_F\) approximation on the same restricted \(3\times3\) candidate coefficient.

**No source production modules are modified; no new physical UV or EFT matching is claimed.** This package is separate from the latest Workstream A changes. The input is the completed `t3_bc_running_one_scale_triangle_v1` report for each model.

## Exact mathematical calculation

With \(\Phi=(S_1,\epsilon S_2^*)\), \(A=\mathcal M_\Phi^2=A^\dagger>0\), and \(U^\dagger A U=\mathrm{diag}(m_X^2,m_Y^2)\), the original-basis Euclidean propagator is

\[
G(t)=(t\mathbf1+A)^{-1},\qquad t=p_E^2.
\]

At zero external momentum the **mass part** of the restricted one-loop triangle can be packaged into a four-index pair kernel:

\[
K_{(ab),(cd)}(M_F)=-M_F^2\int_0^\infty\!dt\,\frac{t}{t+M_F^2}
G_{ac}(t)G_{bd}(t).
\]

The independently evaluated spectral result is

\[
K=(U\otimes U)\operatorname{diag}\big[\Phi(M_F^2,m_i^2,m_j^2)\big]
(U\otimes U)^\dagger,
\]

\[
\Phi(A,B,C)=-A\int_0^\infty\frac{t\,dt}{(t+A)(t+B)(t+C)}.
\]

The spectral \(\Phi\) is evaluated using 55-digit divided differences. The original-basis \(K\) is computed independently using `scipy.integrate.quad_vec` over logarithmic Euclidean momentum. The matrix-valued integral is not reconstructed by inserting the known scalar eigenvalues into the original-basis route.

A second original-basis evaluation expands the propagator using an arbitrary number \(N\) of insertions of the off-diagonal mass \(V=A-\mathrm{diag}(A)\):

\[
G_N(t)=\sum_{n=0}^{N}[-D_0(t)V]^nD_0(t),\quad D_0=(t\mathbf 1+\mathrm{diag}(A))^{-1}.
\]

This checks what is gained by resumming the allowed quadratic instead of deleting it. The calculation is explicitly a **scalar propagator identity**, not a proof of the UV spinor factor or absolute physical neutrino-mass sign. The raw `K` omits \(1/(16\pi^2)\), Yukawas and external vertices.

The Matchete-original seven native quartic tensors are imported via `scripts/AuditT3BCNativeQuarticMassBasis.py`, not replaced by random tensors. Each is rotated at the new **running** \(U\) using `rotate_tensor(...)` and `verify_polynomial(...)`, with a separate Hermitian-conjugate check; the two `110101` contractions must retain rank two. The topology quartic is still treated as a separate baseline source.

## Existing matched coefficient and hard/soft tests

The runner also imports the **existing** Workstream C v22 matched F/Y/X source-polynomial weights \(W_{XX},W_{YY},W_{XY}\). The exact source-matched coefficient is reconstructed from

\[
A_j(M_F)=\frac{1}{16\pi^2 M_F}\sum_{a=XX,YY,XY}W_{a,j}\,\Phi_a(M_F),
\]

using the same evolved heavy-flavour Yukawas, complex quartics, scalar rotation and three Takagi masses as the uploaded v1 result. The exact recomputed \(3\times3\) matrix must reproduce the stored full-mass candidate to relative tolerance \(2\times10^{-8}\). The alternate strict expansion uses

\[
\Phi_{\rm strict}(A,B,C)=\frac{B\ln(B/A)-C\ln(C/A)}{B-C}.
\]

At three distinct scales \(\mu=M_F,M_Y,M_X\), the code independently enforces the v19 restricted hard/soft identity

\[
H(\mu)+S(\mu)=\Phi_{\rm strict},\quad
H=1+\ln(\mu^2/M_F^2),\quad
S=-\left[1-\frac{B\ln(B/\mu^2)-C\ln(C/\mu^2)}{B-C}\right].
\]

This does **not** supply the uncomputed finite parameter, wave-function, EOM, or scalar threshold matching. The v19 identity is restricted to the leading \(1/M_F\) scalar-bubble sector.

The batch also tests an explicitly **synthetic** complex Majorana mass with diagonal phases \((0.19,-0.13,0.27)\) rad. Reusing the same source direct 480-term coefficient, it compares two candidate heavy-spinor rotation prescriptions (`y1@V*` versus `y1@V`). The original all-real benchmark gives exactly zero difference and cannot test the phase convention; the synthetic experiment is a diagnostic only, not a physical spectrum or a replacement for Workstream A's spinor normalization certificate.

## Installation and PowerShell

Install the full ZIP at the repository root:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\t3_bc_mixed_quadratic_kernel_v2.zip" -DestinationPath . -Force
& .\scripts\RunT3BCMixedQuadraticKernelAudit.ps1
```

The runner executes all Python tests and, if successful, performs the full **B and C** source-backed calculation. Requirements: saved v1 running-input reports; `output/general_t3_stage22/matchete_tensors.json`; previous `scripts/AuditT3BCNativeQuarticMassBasis.py`; the Matchete v22/v16/v17/v18 source scripts and their archived source evidence; NumPy, SciPy, mpmath, SymPy. These are all dependencies of existing project calculations. Python source files, upstream RGBeta and Wolfram files are **read-only**.

For a kernel-only diagnostic without archived native tensor or matched-source files (not a full batch completion):

```powershell
python -u -X utf8 -m scripts.RunT3BCMixedQuadraticKernelAudit --input-dir output\studies\t3_bc_running_one_scale_triangle_v1 --skip-native --skip-source-weights --skip-spinor-stress --output-id t3_bc_kernel_partial_debug
```

This produces the distinct status `CompletedSubsetByExplicitSkip`; it must not be described as full source-level validation.

Default outputs (new, never overwritten):

- `output/studies/t3_bc_mixed_quadratic_kernel_v2/B_m1_mixed_kernel.json`
- `output/studies/t3_bc_mixed_quadratic_kernel_v2/C_m1_mixed_kernel.json`
- `output/studies/t3_bc_mixed_quadratic_kernel_v2/summary.json`

To repeat, choose a different ID using `--output-id`. An existing output folder is rejected.

## Verified here against the two uploaded reports

The numerical kernel and hard/soft calculations below were independently executed on `B_m1_running_one_scale.json` and `C_m1_running_one_scale.json` uploaded in this conversation, **without** relying on the local Matchete symbolic source inventory. Therefore the native tensor transport, full matched-source \(3\times3\) strict comparison and complex-spinor stress run require your Windows checkout.

| Quantity | B(-1) | C(-1) |
|---|---:|---:|
| Original-basis integration versus spectral kernel, worst over three MF | 2.34e-13 | 2.27e-13 |
| Error from omitting quadratic-mixing insertions, worst over three MF | 1.065e-3 | 1.032e-3 |
| Error after two quadratic-mixing insertions, representative lowest MF | 2.63e-8 | 2.54e-8 |
| Worst exact versus strict \(\Phi\) mass-function relative correction, all MF and XX/XY/YY | 2.883e-4 | 2.253e-4 |

These are **kernel** relative differences. The actual flavour-assembled \(C_5\) comparison can be larger or smaller due to cancellations. The source-backed strict-\(1/M_F\) \(3\times3\) calculation in this package evaluates that effect directly when run locally.

## Remaining physical gates

1. The full source Matchete UV interaction assembly with independent nonzero mixed quadratic and all baseline/Yukawa native tensors is not regenerated here; a mathematical propagator resummation does not certify the original archived UV matching source in that more general theory.
2. The Matchete ↔ RGBeta physical action dictionary, and heavy-spinor chiral/charge-conjugation rotation, remain dependent on Workstream A.
3. Fully renormalized EFT1/EFT2 finite field/parameter/EOM matching, complete operator closure and threshold-specific RG transport remain unfinished. A positive result **must not** set `physical_C5_certified` to `true`.

The next physically decisive task after this audit is to use Workstream A's certified interaction dictionary to assemble an actual general-mixed-mass Matchete UV Lagrangian and re-run direct and sequential matching, rather than using the archived formula with new parameters as a substitute for that calculation.
