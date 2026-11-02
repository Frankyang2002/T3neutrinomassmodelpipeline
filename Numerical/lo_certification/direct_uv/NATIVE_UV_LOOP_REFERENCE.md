# Workstream 4 — Independent full-mass UV spectral reference (batch 6)

## Scope and certification

This is an **additive** cross-check of the frozen restricted benchmark. It does **not** change the original core, frozen values, Yukawa tensors, effective theories, threshold scales, or RG equations. It is not a cross-package normalisation certificate. The output remains `physical_source_normalisation_certified: false` for all 16 models until the actual native Matchete spinor/CG coefficient map and scalar portal identification have been independently verified.

The previous installed batch, `Numerical/lo_certification/direct_uv/native_amplitude_contract.py`, is a dependency. It was verified locally by the user with 7/7 passing tests and 16/16 frozen source-structure checks. The new batch re-runs those source guards and checks the frozen tensors before performing loop calculations.

## Full-mass independent amplitude

The numerical source action is

\[
\mathcal L \supset -\frac12\psi^TM\psi - \nu_a\psi_i(Y_{ais}\phi_s+Z_{ais}\phi_s^*)+\mathrm{h.c.},\quad
V_2=\phi^\dagger A\phi+\mathrm{Re}(\phi^T B\phi).
\]

With \(\phi_s=(x_s+i y_s)/\sqrt2\), the Cartesian Hessian is

\[
H_R=\begin{pmatrix}
\Re A+\Re B&-\Im A-\Im B\\
\Im A-\Im B&\Re A-\Re B
\end{pmatrix}, \qquad O^TH_RO=\mathrm{diag}(m_{Sr}^2).
\]

Writing \(C_{sr}=(O_{s r}+iO_{n_s+s,r})/\sqrt2\), the vertex into the real scalar \(R_r\) is

\[
G^{(v)}_{air}=Y^{(v)}_{ais}C_{sr}+Z^{(v)}_{ais}C_{sr}^*.
\]

The exact zero-momentum **full-mass** chirality-flipping self-energy is reconstructed without Takagi diagonalisation of \(M\):

\[
\boxed{\quad
M_\nu=\frac1{16\pi^2}\sum_r G_r\, M^\dagger f_r(MM^\dagger)\,G_r^T,
\qquad
f_r(z)=B_0^{\overline{\rm MS}}(0;z,m_{Sr}^2;\mu^2).
\quad}
\]

This follows from the exact identity \(M^\dagger f(MM^\dagger)=U\,\mathrm{diag}(m_k f(m_k^2))\,U^T\) for Takagi \(U^TMU=\mathrm{diag}(m_k)\), but evaluates the left-hand side through a **Hermitian eigendecomposition of \(MM^\dagger\)**. The production engine instead constructs Takagi eigenstates. These are independent numerical methods for the same chosen input action.

The finite function for positive squared masses is

\[
B_0(0;a,b;\mu^2)=1-\ln(a/\mu^2)-g(b/a),\qquad
 g(x)=\frac{x\ln x}{x-1},\quad g(1)=1.
\]

The code first verifies the complete coefficient of the UV pole:

\[
P=\frac1{16\pi^2}\sum_r G_r M^\dagger G_r^T=0.
\]

**Only after this check** it subtracts a common reference scalar mass from \(B_0\), thereby eliminating the cancelling scalar-independent terms *before* floating-point summation. It does **not** discard a nonzero pole or choose a counterterm. The reference squared mass is chosen as the geometric mean of the smallest and largest scalar squared masses; a valid zero-pole amplitude is independent of that arbitrary reference. The amplitude retains exact scalar and heavy masses; there is no large-mass expansion.

All 16 model outputs contain: full complex \(3\times3\) matrix (GeV), three Takagi masses as the singular values (eV), vertex-pair matrices, pole norms, cancellation factor, fermion and scalar spectra, and optional comparison to `Numerical/direct_one_loop/core.py`. A field-rephasing check is available. The comparison to `core.py` is conditional and may reveal large relative errors when individual terms nearly cancel; *inspect the absolute and natural-amplitude errors as well as the matrix-relative error.* No arbitrary agreement tolerance is used to silently recertify a case.

## Source-level scalar portal issue — specifically verified, not silently fixed

At the currently inspected source revision (3 November 2026), `Lagrangian/interactions/ScalarPotential.wl` constructs portal couplings dynamically as `lambdaH1InvN` and `lambdaH2InvN`, whereas `Lagrangian/model/T3Fields.wl` separately defines unindexed `lambdaH1`, `lambdaH2`. `Lagrangian/interactions/T3Topology.wl` assembles `BuildGeneralScalarPotentialCandidates`. Meanwhile `scripts/BuildT3DirectOneLoopInputs.py` reads the frozen `lambdaH1`, `lambdaH2` to construct the provisional scalar matrix. An equality between these names, or a factor/sign conversion, is **not proved** just by these source structures. The batch records this gap, never sets the two equal, and does not modify numerical input matrices.

To obtain an actual physical portal dictionary, isolate one registered native `lambdaHiInvN` operator, substitute the physically normalised Higgs neutral component, and compare the resulting Cartesian scalar Hessian with the frozen parameter definition. In non-singlet representations, include every independent invariant and the native CG normalisation. The native operator-export script below provides an inspectable `LInt` inventory for the two representative classes, but **does not yet perform** that differential Hessian measurement.

## Optional live Matchete native evidence

`scripts/RunLODirectUVNativeOperatorExport.wl` builds the original restricted `A_+0` and `B_-1` Matchete UV Lagrangians using the existing `BuildT3Lagrangian`, verifies `FullValidation`, and exports the full symbolic assembled `LInt` plus presence profiles of portal and Yukawa function names. It runs neither `Match` nor RGEs. It refuses an existing output and requires all experimental extension environment flags not be `1`.

If supplied via `--native-export`, the Python runner validates the exported **actual source hashes** against the currently installed Lagrangian. It rejects stale evidence and a missing pair of topology operators. Even a passing native inventory does **not** fix the concrete `GammaCC` phase, RGBeta epsilon/triplet basis, the full y1/y2/MF transport, or any signed portal coupling mapping. These remain outstanding.

## Installation and execution (Windows PowerShell)

Extract the single ZIP over the **root** of `C:\thesistemp\T3neutrinomassmodelpipeline` (do not extract into a nested folder). The ZIP contains complete new files and does not replace the previous batch.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

.\.venv\Scripts\python.exe -m pytest `
    tests\test_lo_direct_uv_native_loop_reference.py -q

.\.venv\Scripts\python.exe `
    scripts\RunLODirectUVNativeLoopReference.py `
    --project-root . `
    --rephase `
    --output-dir output\general_observable\lo_certification_direct_uv\native_uv_reference_01
```

The default command reads `output/general_observable/lo_direct_one_loop/provisional_inputs/<SIGNED_MODEL>.json`, verifies the existing source-operator contract and all 16 input shapes, recomputes the UV mass **independently**, then compares it to the current production direct-UV one-loop engine (read-only). It writes:

- `output/general_observable/lo_certification_direct_uv/native_uv_reference_01/all_16_native_uv_reference_status.json`
- `output/general_observable/lo_certification_direct_uv/native_uv_reference_01/all_16_native_uv_reference_summary.csv`
- `output/general_observable/lo_certification_direct_uv/native_uv_reference_01/models/<SIGNED_MODEL>.json` (16 files).

The path must be new or empty. If rerunning, use `native_uv_reference_02` etc. The `--no-production-comparison` option runs only the independent reference. The `--rephase` option repeats each calculation after a complex field-basis change and reports both natural-scale and matrix-relative errors.

**Optional** (only if Matchete can be loaded by `wolframscript` locally):

```powershell
wolframscript -file scripts\RunLODirectUVNativeOperatorExport.wl `
    output\general_observable\lo_certification_direct_uv\native_operator_export_01\native_operators.json

.\.venv\Scripts\python.exe `
    scripts\RunLODirectUVNativeLoopReference.py `
    --project-root . `
    --native-export output\general_observable\lo_certification_direct_uv\native_operator_export_01\native_operators.json `
    --output-dir output\general_observable\lo_certification_direct_uv\native_uv_reference_02
```

If the Wolfram execution is blocked or times out, retain the independent calculation and classify the native operator measurement as missing. Do not replace it with synthetic positive evidence.

## Validation limits

Eleven local synthetic numerical tests have been executed successfully (including an all-16 synthetic model batch). The **actual 16 frozen numerical inputs and native Matchete** were not available in this environment: the local commands above must establish those results. The existing source snapshots were inspected from the read-only GitHub repository, not modified there.

### Workstream boundaries

At the frozen input scale \(\mu_0=10^7\) GeV this reference is a full UV one-loop calculation **conditioned on the supplied component action**. It is not an EFT1 boundary at \(\mu_F=10^5\) GeV. Keep the sequence UV \(\to\) integrate out **all** heavy fermions \(\to\) EFT1 with **both** heavy scalars and SM fermions \(\to\) integrate out **both** scalars \(\to\) SMEFT. The all-heavy-fermion Gaussian LLSS bridge is a separate tree-level consistency object; hard one-loop fermion matching, scalar matching and complete UV/EFT/SMEFT agreement are excluded.
