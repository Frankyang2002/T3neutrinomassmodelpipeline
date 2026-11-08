# A(-4): five-monomial RGBeta diagnostic and conditional ODE comparison

## Scope

Checked against GitHub `main` revision `54857efafc7261f0e4c2bd0802f23e61dc3c165a`, 8 November 2026. This batch introduces **new diagnostic files only**; it does not change the native RGBeta implementation, published-source Wilson anomalous dimension, Matchete input, benchmark definitions, or SMEFT matching.

There are three different claims, and they must not be conflated:

1. **Verified from the saved 55-probe scalar-only numerical audit:** exactly five independent native-minus-reference quadratic monomials reproduce the 34 failing unit/pair directions. The scripts **recompute** this rather than trusting status strings.
2. **Conditional diagnostic only:** removing those five scalar-only differences from a *copy* of the native beta allows a controlled numerical comparison of renormalisable A(-4) EFT1 running in the *declared* 12-real-field reference basis. It does **not** prove native RGBeta wrong; the Matchete-to-RGBeta invariant map and full gauge/Yukawa quartic beta are not independently certified.
3. **Still unresolved:** installed Matchete's signed `GammaCC` conversion to the chosen master Weyl bilinear. Both `rho=+1` and `rho=-1` source reparametrisations remain viable without an independent signed connected amplitude or explicit native convention.

## Exact native-minus-reference scalar-loop polynomial

Define `delta_beta = beta_native - beta_reference`, with every beta in `16*pi^2*d(lambda)/d ln(mu)` convention. With `a=lambdaS2Adj`, `b=lambdaH2Adj`, and `z=lambdaT3`, the only nonzero components are:

- `delta_beta_lambdaS2 = 8*a^2`
- `delta_beta_lambdaS2Adj = 14*a^2 + 2*b^2`
- `delta_beta_lambdaH = 2*b^2`
- `delta_beta_lambdaH2 = 4*b^2`
- `delta_beta_lambdaH2Adj = 10*a*b - (16/3)*|z|^2`

The rest of the ten scalar-only beta-coordinate discrepancies vanish within the previous numerical audit precision. The diagnostic alternative is `beta_diagnostic = beta_native - delta_beta` in **only** these five components, leaving all other 12 of the 17 EFT1 beta expressions unchanged. It uses `Bar[lambdaT3]*lambdaT3` in the existing safe InputForm parser for complex conjugation. This is a deliberately nonphysical reference-comparison payload, **not a new production beta**.

## Outputs

### `scripts/RunAminus4QuarticFiveMonomialBatch.py`

Reads `output/a_minus4_eft1_batch2/A_m4_projected_raw.json`, `A_m4_projection_audit.json`, `A_m4_eft1_tensor_basis.json`, plus the three existing signed-Weyl witnesses.

Writes (into the chosen new output folder):

- `five_monomial_report.json`: independently evaluated 55 unit/pair probes, five coefficients, five random mixed-quartic checks, complex phase check, SHA-256 provenance, simple fixed-background sensitivity, and signed-spinor status.
- `native_global_projection_compatible.json`: native projected 17-beta baseline, only with the metadata needed by the existing ODE evaluator.
- `scalar_only_reference_DIAGNOSTIC_beta.json`: a **separate nonphysical** alternative passed through 55 scalar-only reference tests and five random complex-point tests. Never substitute it silently as a physical result.

### `scripts/RunAminus4ConditionalQuarticODE.py`

Uses the same archived initial state at `mu_F` and integrates the **existing** 17-coupling ODE with each beta source, with both `S1` and `S2` active down to `mu_S`. It saves full 45-Weyl/12-real-scalar running tensors for **both** and compares the ten quartic coordinates:

- `conditional_ode_summary.json`
- `coupling_differences.csv`
- `native/A_m4_eft1_tensors.npz`, `native/A_m4_eft1_tensor_basis.json`
- `reference_diagnostic/A_m4_eft1_tensors.npz`, `reference_diagnostic/A_m4_eft1_tensor_basis.json`

These are **conditional renormalisable EFT1 trajectories**, not physical neutrino masses or SMEFT Wilson coefficients.

## PowerShell — run from project root

```powershell
Set-Location "C:\thesistemp\T3neutrinomassmodelpipeline"

Expand-Archive `
  -Path "C:\Users\popok\Downloads\Aminus4_QuarticFiveMonomial_ConditionalODE_Batch.zip" `
  -DestinationPath "C:\thesistemp\T3neutrinomassmodelpipeline" `
  -Force

python -m unittest discover -s tests -p "test_a_minus4_quartic_five_monomial_batch.py" -v

python scripts\RunAminus4QuarticFiveMonomialBatch.py `
  --output-root "output\a_minus4_next_batch\quartic_five_monomial"

python scripts\RunAminus4ConditionalQuarticODE.py `
  --input-root "output\a_minus4_next_batch\quartic_five_monomial" `
  --output-root "output\a_minus4_next_batch\conditional_quartic_ode"
```

The commands intentionally **refuse to overwrite** existing diagnostic JSONs and ODE outputs. Use a new output folder to rerun.

## Next physical gates

- Establish the installed RGBeta quartic registry's exact native *component-level* contraction and its map to the real fourth-derivative tensors; do not infer a correction from the five-monomial discrepancy alone.
- Derive the installed Matchete `GammaCC` signed Weyl amplitude or an equivalent fully source-fixed connected four-point UV amplitude. The existing abstract relations only fix it up to sign and field/coupling redefinitions.
- Continue with the required `UV -> integrate out all F -> EFT1(SM+S1+S2) -> integrate out all scalars simultaneously -> SMEFT` order. Neither diagnostic output is a physical SMEFT C5.
