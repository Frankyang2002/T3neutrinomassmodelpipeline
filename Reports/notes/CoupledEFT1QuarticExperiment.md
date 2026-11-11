# Coupled EFT1 quartic-only LLSS–LLHH test

**Status:** Exploratory. Intended to test the numerical coupling without
altering the authoritative `corrected_llss16_v4` run.

## Files in this ZIP

- `scripts/RunT3CoupledEFT1Quartic.py` — complete new executable module
- `Reports/notes/CoupledEFT1QuarticExperiment.md` — these instructions

## Physics scope

The script reuses the existing 16-model fermion-tree source, endpoint HHSS
quartic tensor, finite joint-scalar bubble, renormalisable UV/EFT1 running,
and SMEFT running. It replaces **only the EFT1 LLSS leading-log evaluation**
with a coupled numerical LLSS↔LLHH quartic-sector ODE.

Implemented:

    16*pi**2 d C_LLHH / d ln(mu) = contract(lambda_HHSS, C_LLSS)
    16*pi**2 d C_LLSS / d ln(mu) = contract(lambda_HHSS, C_LLHH)

C_LLHH(muF) = 0 *within this retained LLSS approximation*. Both tensors
retain all flavour, weak and scalar-component indices throughout integration.
Quartics are fixed to the existing EFT1 endpoint; the original renormalisable
trajectory is **not** integrated simultaneously with Wilson coefficients.

Not implemented: other five master-anomalous-dimension interactions,
gauge/Yukawa Wilson corrections, HHHH/SSSS quartic feedback,
complete operator inventory, hard one-loop fermion matching, convention
certification. The `W=1` result is thus *NOT the full EFT1 numerical solution*.

The existing scalar finite bubble is deliberately kept identical to the
fixed-order baseline. The only comparison is fixed-order versus the
coupled quartic-only Wilson evolution.

## Install (PowerShell)

Run at the repository root. Make sure ZIP is in Downloads:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -Path "$env:USERPROFILE\Downloads\T3_Coupled_EFT1_Experiment.zip" -DestinationPath . -Force
```

## First tests (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

# Independent mathematical solver test: no project inputs needed
python -X utf8 -m scripts.RunT3CoupledEFT1Quartic --self-test

# Confirm all existing fitted input files are present
python -X utf8 -m scripts.RunT3CoupledEFT1Quartic --source-run thesis16_merged_v1 --run-id coupled_A_m4_v1 --models A_-4 --check-only

# Actually compute the A(-4) experimental cases (16 U/E/W/S combinations)
python -X utf8 -m scripts.RunT3CoupledEFT1Quartic --source-run thesis16_merged_v1 --run-id coupled_A_m4_v1 --models A_-4
```

Do not rerun the experiment with the same run ID. If needed choose a fresh ID.

## Output location

`output/studies/general_scalar_extended/experimental_coupled_eft1_quartic/coupled_A_m4_v1/`

Inspect the run's `summary.json` and model's `cases.json`, especially
`1111` and `1101`. In `cases.json` each U/E endpoint includes:

- `c5_induced_wilson_log` — integrated quartic-sector response
- `c5_fixed_order_llss_log` — original fixed-order leading log
- `coupled_log_over_fixed_order_log_rel` — relative matrix difference
- `c_llss_endpoint_relative_drift` — running-induced LLSS change
- `first_order_reference_agreement_rel` — first-order parity check
- `solver_coupled_nfev` — number of ODE function evaluations
- `coupled_weinberg_projection_residual` — diagnostic of non-Weinberg pieces

The script fails deliberately if the independently computed fixed-order
derivative does not match the existing LLSS leading-log calculation.

## Extending to all 16 later

Remove `--models A_-4` and use a new `--run-id`. This should only be
done after checking the one-model output.

## Expected difference

The generated Weinberg coefficient is already a one-loop result. Its first
correction due to LLHH↔LLSS back-reaction appears at third order in
the quartic-mixing block; LLSS itself changes at two-loop order.
Therefore a small difference from the fixed-order leading log is expected.

A numerical solver alone does not supply missing one-loop matching or other
anomalous-dimension terms. This experiment establishes integrator behaviour,
not complete physical EFT1 correctness.

## Validation carried out in this chat

- Python syntax compilation: passed
- Isolated n=1 analytic sinh/cosh ODE: passed
- n=8 tensor-shape numerical integration: passed
- Full T3 A(-4) using repository data: **not run in this environment**
  (the repository/data are on the user's local machine).
- GitHub: **read only**, no modifications.
