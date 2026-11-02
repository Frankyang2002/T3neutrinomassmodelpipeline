# Workstream 4 — native `lambdaT3` scalar mixing (restricted UV)

## Scope

Read-only, additive, conditional physical dictionary. The native Wolfram script
builds the **original restricted** Matchete Lagrangian for representatives
`A_+0`, `B_-1`, `C_-1`, `D_+0`, `E_+0`. A second script checks every one of the
**16** unchanged frozen benchmark scalar mixing matrices.

**This does not establish a complete physical source-to-Weyl normalisation.**
No benchmark parameters, one-loop matrices, RGE data or EFT kernels are
modified. The outputs are for audit only.

## Source conventions checked

The native `BuildT3MixingCandidates` term is

\[
\mathcal L_{\rm int}\supset
 +\lambda_{T3}\,C_{ij st}\,H_iH_jS_{1s}S_{2t}^{*}
 +\mathrm{h.c.}
\]

Its native registered tensor `T3MixCG` has flags
`(H,H,S1,Bar[S2]) = (False,False,False,True)`, symmetric Higgs indices, and
**no explicit `1/2`** on the non-Hermitian-conjugate topology operator.
The exporter requires **both** non-H.c. and H.c. expressions to appear in the
*assembled*, validated `LInt` and reads the real registered tensor using
`Matchete`GetCGTensor[T3MixCG]`, never a synthetic CG.

Taking \(H^0=v/\sqrt{2}\) and \(V=-\mathcal L\) for static fields gives the
restricted mixing potential

\[
V_{\rm mix} =
-\frac{v^2}{2}\lambda_{T3}\sum_{st} C_{00st} S_{1s} S_{2t}^{*}
+\mathrm{h.c.}
\]

and, with \(V_2=\phi^\dagger A\phi\),

\[
A_{d_1+t,\;s}= -\frac{v^2}{2}\lambda_{T3} C_{00st},
\qquad A_{s,\;d_1+t}=A_{d_1+t,\;s}^{*}.
\]

The indices `00` mean the two *neutral weight-basis* Higgs components;
Python index `[1,1]` (0-based) follows the existing input builder. It is
**not** assumed to be `[1,1]` before rotating the raw Matchete axes.

This derivation does **not** assign the same physical numerical value to
Matchete's `lambdaT3` and RGBeta's `lambdaT3`: that cross-package mapping
remains a separate native measurement.

## Corrected audit file lookup and provenance

The original Python audit incorrectly required three files under
`Reports/output/full/multiplicity/`. The current repository stores the files
under the **committed canonical** directory
`output/general_observable/lo_direct_one_loop/source_audits/`:

- `component_support_audit.json`
- `matchete_generator_analysis.json`
- `cg_covariance_audit.json`

The frozen component inputs also record this canonical path (Windows
backslashes are normalised when checked) and SHA256 digests of the first
two audits. The fixed audit requires the source files to match each frozen
model's stored provenance and also checks the benchmark SHA256. It never
regenerates or silently substitutes CGs. If these requirements fail, the
audit blocks rather than assuming the native mixing comparison passed.

The original five-class Matchete export is reusable; there is **no need**
to rerun `wolframscript` after this Python-only fix.

## Independent checks

1. Live native `GetCGTensor` for 5 classes, from a successful original
   `BuildT3Lagrangian` and `CheckLagrangian`. Strict source digests prevent
   stale-output reuse. The four opt-in flags are forced to `0` in the child
   Wolfram process, with inherited values recorded.
2. Every full native `(H,H,S1,Bar[S2])` tensor is compared against the
   historical source component-support audit in the **original registered
   Matchete basis**. There is no allowance for a free overall sign or phase.
3. The CG is rotated with the recorded unitary eigenvector matrices for
   dimensions 1, 2 and 3. The `Bar[S2]` index uses conjugate rotation.
4. Nonzero neutral Higgs mixing entries must connect scalar fields of equal
   electromagnetic charge, for all 16 alpha assignments.
5. Each complex off-diagonal block of the *frozen* Hermitian scalar matrix
   `A` is compared against the source-derived expression above, using the
   frozen `lambdaT3` and VEV. The restricted holomorphic `B` must be zero.
6. An independent finite-difference real Cartesian Hessian is reconstructed
   for the complex scalar action, without calling the numerical loop core.
7. Tests include conjugate/sign errors, missing native operators, stale
   source digests, invalid SU(2) rotations, non-Hermitian frozen inputs and
   an end-to-end **synthetic** sixteen-model fixture.

No heavy-mass expansion or further neutrino-mass calculation is undertaken.

## Project-relative paths (five complete additive files)

- `scripts/RunLODirectUVNativeMixingTensorExport.wl`
- `Numerical/lo_certification/direct_uv/native_mixing_dictionary.py`
- `scripts/RunLODirectUVMixingDictionary.py`
- `tests/test_lo_direct_uv_native_mixing_dictionary.py`
- `Numerical/lo_certification/direct_uv/NATIVE_MIXING_DICTIONARY.md`

## Commands (PowerShell, from repository root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

.\.venv\Scripts\python.exe -m pytest `
  tests\test_lo_direct_uv_native_mixing_dictionary.py -q

wolframscript -file `
  scripts\RunLODirectUVNativeMixingTensorExport.wl `
  output\general_observable\lo_certification_direct_uv\native_mixing_tensor_01\native_five_class_mixing_tensors.json

.\.venv\Scripts\python.exe `
  scripts\RunLODirectUVMixingDictionary.py `
  --project-root . `
  --native-json output\general_observable\lo_certification_direct_uv\native_mixing_tensor_01\native_five_class_mixing_tensors.json `
  --output-dir output\general_observable\lo_certification_direct_uv\native_mixing_dictionary_02
```

The supplied native JSON is reused; the Python command writes only
`native_mixing_scalar_A_dictionary.json` in the specified new/empty folder.
Both refuse to overwrite existing outputs. If the native export fails, do
not run the dictionary solver on missing or stale JSON.

## Certification boundaries and next step

The native mixing CG and conditional scalar-Hessian map can become
source-verified after running both native and Python scripts **without**
certifying the physical Matchete/RGBeta relative coupling normalisation.
The remaining relevant native measurement is the RGBeta `delS2` mixing
invariant expressed in the **same** concrete SU(2) axes as Matchete's CG,
including the signed coupling/Feynman-rule prefactor. Then establish the
common signed `GammaCC`, `y1`, `y2`, and `MF` Weyl mapping.

EFT order stays unchanged:

UV (\(\mu_0=10^7\,\mathrm{GeV}\)) \(\to\) integrate out **all BSM
fermions** at \(\mu_F=10^5\,\mathrm{GeV}\), retaining both BSM scalars and
all SM fermions \(\to\) EFT1 \(\to\) integrate out **both** BSM scalars at
\(\mu_S=10^3\,\mathrm{GeV}\) \(\to\) SMEFT.
