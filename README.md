# Candidate scalar endpoint to one-scale neutrino mass: B(-1)/C(-1)

## New paths (DO NOT overwrite root README.md)

- `scripts/RunT3QuadraticEndpointC5.py`
- `tests/test_t3_quadratic_endpoint_c5.py`

All older scripts and outputs are unchanged. Extract ZIP into the repository root.

## Prerequisites

Existing Stage 33 boundary, frozen fitted-reference files, experimental Matchete mass-basis certificates and **saved** Stage quadratic EFT1 flow from `scripts.RunGeneralT3QuadraticEFT1`. This runner deliberately does not launch Wolfram or refit Yukawas.

For missing endpoint flow outputs, generate the two candidates first (Wolfram required):

```powershell
python -m scripts.RunGeneralT3QuadraticEFT1 --model B_-1 --mu2-fraction 0.01
python -m scripts.RunGeneralT3QuadraticEFT1 --model C_-1 --mu2-fraction 0.01
```

Use already saved outputs when available; do not regenerate unnecessarily.

## Run

```powershell
python -X utf8 -m pytest -q tests/test_t3_quadratic_endpoint_c5.py
python -m scripts.RunT3QuadraticEndpointC5 --all
```

Output: `output/studies/general_scalar_extended/quadratic_endpoint_c5/summary.json` and separate B/C reports. For the zero-extra-quartic control:

```powershell
python -m scripts.RunT3QuadraticEndpointC5 --all --zero-extra-quartics
```

**Note:** Running `--zero-extra-quartics` into the default output directory replaces this *new experimental output folder's* B/C files; use `--output-dir output/studies/general_scalar_extended/quadratic_endpoint_c5_zero_extra` to preserve both.

## What the calculation actually does

1. Reads saved candidate running scalar masses, extra mixed quadratic coefficient and quartics at the EFT1 endpoint.
2. Constructs and verifies the B/C complex-doublet mass eigenbasis.
3. Loads SHA256-checked frozen fitted-reference Yukawas/fermion masses (not Yukawa running on the candidate background).
4. Evaluates the already extracted **direct one-scale** Matchete ordered coefficient using the endpoint eigenmasses and mixing matrix.
5. Compares the original `lambdaT3`-only C5 with the two extra C5-active quartics, mapped by the existing source-level sign convention `lambda_matchete = -lambda_rgbeta`.
6. Records Takagi masses and the relative Frobenius norm of the response.

**Not implemented/certified here:** full sequential EFT1 matching or Wilson running, complete field-component/Matchete-RGBeta normalisation, genuine physical predictions, self-bilinear A/D/E threshold treatment. Reports explicitly say so. The nonzero-extra-quartic result is a *conditional sensitivity probe*.
