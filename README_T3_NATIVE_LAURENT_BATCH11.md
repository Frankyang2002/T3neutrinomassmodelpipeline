# T3 Batch 11 — Native F-hard + joint-scalar Laurent closure

A(−4) and D(−2) only. This is an **additive local-only implementation** following successful Batch 9 and Batch 10. It does not edit the existing stages or GitHub and does not invent physical matching terms.

## Implemented

1. **Read the actual saved full-general Matchete scalar-soft expressions**, not a reconstructed bubble proxy:
   - A(−4): `output/parallel_t3/full_general_final_one_batch/thesis_a4_fresh_f_v3/A_-4/A_scalar_soft_only_unrenormalized.txt` (10 native soft Weinberg terms).
   - D(−2): `output/parallel_t3/full_general_final_one_batch/thesis_physical_scalar_continuation_v1/D_-2/D_scalar_soft_only_unrenormalized.txt` (68 native terms in the physical split-real `T3Scalar2Minus`/`T3Scalar2Plus` basis).
   - Verify the archived native success status and model before using either expression.
2. **Extract the exact Laurent pole and regular parts** with SymPy, explicitly recognizing only the source's `T3EFT1W12`, `T3EFT1WExtra`, `lambdaT3`, `lambdaT3ExtraQ020101I1`, heavy mass, `hbar`, `ϵ`, and `μbar2` structures. Other operators and couplings cause a hard error; none are silently set to zero.
3. Insert the existing three-flavour F-threshold flavour contractions, using the Batch 1 F-endpoint Yukawas and positive heavy masses:

   \[
   W_{12}^{pq}=\sum_r\frac{y_{1,pr}y_{2,qr}}{M_r},\qquad
   W_{\mathrm{extra}}^{pq}=\sum_r\frac{y_{1,pr}(y_{2\mathrm{Conj}})_{qr}}{M_r}.
   \]

   D's benchmark has `y2Conj=0`, which is explicitly checked; it does **not** drop its nonzero extra scalar coupling. For the native scalar-soft expression, the numerical identification `lambdaT3ExtraQ020101I1 = lambdaT3GeneralQ020101I1` is a **conditional convention assumption**, not an independently verified coupling normalization. Nonzero D extra Yukawa remains unsupported by Batch 10's hard source parser, and will block.
4. **Test native hard/soft pole cancellation** at common frozen one-loop coupling coordinates. For the ordered source amplitudes \(H\) and \(S\):

   \[
   P_H+P_S \stackrel{?}{=}0,
   \qquad \frac{\partial S_{\mathrm{reg}}}{\partial \ln\bar\mu}=2P_S.
   \]

   These are source-level tests, not MS-bar counterterm certification. The code symmetrizes each 3×3 matrix with an ordinary transpose, never a Hermitian conjugate.
5. **Assemble the explicit matching-scale logarithm exactly once**, derived from the native soft logarithm:

   \[
   \Delta C^{\mathrm{source}}_{\mathrm{EFT1}}
      =-2P_S\ln\frac{\mu_S}{\mu_F},\qquad
   C^{\mathrm{source}}_{\mathrm{finite}}
      =H_{\mathrm{reg}}(\mu_F)+S_{\mathrm{reg}}(\mu_S)
       +\Delta C^{\mathrm{source}}_{\mathrm{EFT1}}.
   \]

   It verifies independence under simultaneous independent changes of the reference F/S matching scales with couplings and scalar masses held fixed. It does **not** add Batch 9's scalar log/B0 or Batch 10's conditional extra-HHSS reference tensors a second time.
6. Optionally feed the **explicitly unverified** source-sector finite matrix through the existing SM/Weinberg RGE and NuFIT-Gaussian/PMNS machinery for a numerical sensitivity only: `--illustrative-smeft`. The default run does not fit or claim a physical mass.

## What's **not** finished

This pole-closed source sector is **not a certified physical \(C_5\)**. In particular: the absolute Matchete GammaCC → left-Weyl factor/sign is open; the D native/RGBeta extra-quartic factor is not independently derived; the **full** EFT1 dim-5 tensor ADM and all finite F/UV hard operators have not been included; the native D scalar-source inventory does not contain every additional HHSS channel from the independent reference-basis registry; and a fully source-consistent 18-real-Yukawa refit has not yet been performed. Numerical pole cancellation is necessary but not sufficient for physical one-loop matching.

Strict fixed-order one-loop treatment: freeze the dimensionless scalar couplings at \(\mu_F\) inside both one-loop graphs; use the actual scalar mass spectrum at \(\mu_S\). Differences due to running those couplings inside already-one-loop terms begin at two loops. All three F flavours are removed at \(10^5\) GeV, with both scalar multiplets retained until simultaneous removal at \(10^3\) GeV. Low energy is 100 GeV.

## ZIP contents / project-relative paths

| Path | Full file |
|---|---|
| `Numerical/parallel_t3/complete/two_model_native_laurent_batch11.py` | Source reader, symbolic Laurent extraction, three-flavour substitution, strict-order source log, pole and scale checks |
| `scripts/RunT3NativeLaurentBatch11.py` | Read-only input preflight and two-model runner with new output owner |
| `tests/test_t3_native_laurent_batch11.py` | A analytic pole/log test, matching-scale invariance, D complex flavour and unknown-coupling tests, plus saved native source parse tests |
| `README_T3_NATIVE_LAURENT_BATCH11.md` | This document |

No Batch 9, Batch 10, or existing physics file is replaced.

## PowerShell — install, test, run

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "$HOME\Downloads\t3_native_laurent_batch11.zip" `
  -DestinationPath . -Force

.\.venv\Scripts\python.exe -m pytest -q `
  tests\test_t3_native_laurent_batch11.py

.\.venv\Scripts\python.exe -X utf8 `
  -m scripts.RunT3NativeLaurentBatch11 --all --preflight

.\.venv\Scripts\python.exe -X utf8 `
  -m scripts.RunT3NativeLaurentBatch11 --all --debug
```

Default results:

- `output/parallel_t3/thesis_native_laurent_batch11_v1/A_-4/result.json`
- `output/parallel_t3/thesis_native_laurent_batch11_v1/D_-2/result.json`
- `output/parallel_t3/thesis_native_laurent_batch11_v1/summary.json`

Expected status when source algebra and thresholds agree:
`NativeHardSoftLaurentAndScaleLogClosed_ConditionalOperatorConvention`.

To rerun after a failed attempt, choose a **new output folder**, e.g.:

```powershell
.\.venv\Scripts\python.exe -X utf8 `
  -m scripts.RunT3NativeLaurentBatch11 --all --debug `
  --output output/parallel_t3/thesis_native_laurent_batch11_retry_v2
```

To point to a different previously successful Batch 10 run, supply `--batch10 <relative-path>`; the default is the confirmed `output/parallel_t3/thesis_hard_general_hhss_batch10_fix_v2`.

## Build validation

In the construction environment, four isolated analytic/lexer tests passed and two tests requiring the original A/D native output files skipped. The current saved native source filenames, documented success statuses, and D source coupling-symbol inventory were separately inspected read-only from the repository. Full numerical evaluation against your checked-out Matchete sources still requires the Windows commands above. This distinction is important: the ZIP contains code, **not** a claimed successful run of those native inputs on your machine.
