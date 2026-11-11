# Simultaneous EFT1 RGBeta + seven-term Wilson integration: A(-4)

**Status:** Standalone exploratory A(-4) experiment. No existing source is overwritten. No GitHub changes.

## Project files in this ZIP

- `scripts/RunT3SimultaneousEFT1A4.py` — full new runner
- `Reports/notes/SimultaneousEFT1A4.md` — this note

### Prerequisites

This script reuses existing local modules and the two experimental files already used successfully:

- `scripts/RunT3CoupledEFT1SevenTerm.py` (the **ParserFix** version)
- `scripts/RunT3CoupledEFT1Quartic.py`

Also requires the project's ordinary fitted source campaign `thesis16_merged_v1` and cached EFT1/Matchete input files, as checked by `--check-only`.

## What changes numerically

For the EFT1 stage only, DOP853 integrates in **one** initial-value problem:

\[
 X(t)=(X_{\mathrm{ren}}(t),\operatorname{Re}(C(t))/N,\operatorname{Im}(C(t))/N),\qquad t=\ln\mu,
\]

\[
 \frac{dX_{\mathrm{ren}}}{dt}=\frac{1}{16\pi^2}\beta_{\mathrm{RGBeta}}(X_{\mathrm{ren}}),
 \quad
 \frac{dC}{dt}=\frac{1}{16\pi^2}\beta_C^{(7)}(C,X_{\mathrm{ren}}).
\]

Here \(N=\|C(\mu_F)\|_F\) is only a numerical scaling factor.

- UV running, fermion-tree LLSS matching, joint scalar-pair finite matching, and SMEFT running are reused unchanged.
- Both scalar multiplets remain active throughout EFT1 and are integrated out jointly at the scalar threshold.
- The full *retained* nine-lepton-Weyl-field Wilson tensor is integrated, including initially zero components; all seven master beta terms are evaluated from the instantaneous couplings.
- The quark Yukawa contribution to the Higgs scalar anomalous dimension is retained from the previous seven-term experiment.
- HHSS quartics are built from the original source tensors as a **linear** function of running quartics and checked against the original builder at both endpoints. Non-HHSS HHHH/SSSS quartics use the previous A(-4) explicit baseline-potential ansatz.
- The `E=0` ablation holds the renormalisable state constant while still integrating the Wilson tensor; `E=1` is the simultaneous run. The `W=0` ablation retains the original finite scalar bubble without the induced Wilson mixing contribution.

**Limits:** This is not complete independently certified one-loop neutrino matching. There are no independent hard one-loop fermion threshold terms or dimension-five operator classes outside \(\psi^2\phi^2\). The full 45-Weyl-field operator basis is not integrated. GammaCC-to-Weyl and potential-basis conventions remain provisional. This experiment checks the desired integration methodology, not the correctness of all physical matching inputs.

## Install (PowerShell)

Download the ZIP to Downloads, then:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -Path "$env:USERPROFILE\Downloads\T3_Simultaneous_EFT1_A4.zip" -DestinationPath . -Force
```

## Run (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m scripts.RunT3SimultaneousEFT1A4 --self-test
python -X utf8 -m scripts.RunT3SimultaneousEFT1A4 --source-run thesis16_merged_v1 --run-id simultaneous_A_m4_v1 --models A_-4 --check-only
python -X utf8 -m scripts.RunT3SimultaneousEFT1A4 --source-run thesis16_merged_v1 --run-id simultaneous_A_m4_v1 --models A_-4
```

Use a new `--run-id` if the output directory exists. The runner will not resume or overwrite existing data.

## Result files

```
output/studies/general_scalar_extended/experimental_simultaneous_eft1_seven_term/simultaneous_A_m4_v1/summary.json
output/studies/general_scalar_extended/experimental_simultaneous_eft1_seven_term/simultaneous_A_m4_v1/A_m4/cases.json
```

Inspect the `1111` case (UV/EFT1/Wilson/SMEFT all enabled). Key fields:

- `simultaneous_vs_frozen_7term_rel`: combined vs previous frozen-background master ODE
- `simultaneous_vs_fixed_log_rel`: combined vs original LLSS leading log
- `separate_rgbeta_endpoint_max_block_rel`: independent reference check for renormalisable endpoint
- `initial_scalar_quartic_source_parity_rel`: master scalar beta vs existing LLSS source at muF
- `hhss_linear_basis_parity_rel`: fast linear HHSS reconstruction vs existing builder at muS
- `combined_ode_nfev`: DOP853 RHS evaluations
- `combined_ode_weinberg_projection_residual`: diagnostic for generated tensor not exactly Weinberg-shaped

The runner deliberately blocks if the direct HHSS reconstruction, initial master quartic contraction, or renormalisable trajectory parity checks fail.

Expected: E=0 should reproduce the earlier frozen-coupling seven-term result, whereas E=1 may differ because gauge, Yukawa and quartic couplings vary continuously along the EFT1 interval.

## Validation in this chat

- Python syntax compilation: passed.
- Complex Wilson-tensor packing roundtrip: passed in isolation.
- Analytic toy with an exponentially running renormalisable coupling and coupled Wilson subsystem: passed.
- Full A(-4) fitted source + RGBeta payload run: **not run here**; run locally using the commands above.
