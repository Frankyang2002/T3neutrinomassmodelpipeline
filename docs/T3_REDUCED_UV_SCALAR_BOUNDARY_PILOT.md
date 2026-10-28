# A(-2) UV-running scalar mass-boundary feasibility pilot

## Scope

The completed frozen-UV Sobol/A(-2) comparison achieved 32/32 accepted individual fits. The remaining original running-UV `small_uv_running_5e-9` scenario was screened out for all 16 models. The diagnostic records had scalar masses about 30–169 TeV at the **fixed scalar threshold** of 1 TeV, above the retained 0.5–5 TeV mass window; some models also failed the fermion-first hierarchy diagnostic.

This new script checks **one model, A(-2)** without modifying any earlier campaign files. It uses the source-verified candidate one-loop UV and EFT1 RGBeta equations of the current repository, including all five A(-2) additional quartic sidecar coefficients starting at **exactly zero**, and the existing fermion-threshold projection. No additional quadratic term is added. The topology coupling is fixed at `lambdaT3(1e7 GeV)=5e-9`.

The pilot varies **only** the two *UV scalar squared-mass input parameters* at `1e7 GeV`; all Yukawa matrices, fermion masses, other scalar quartics, thresholds and RGEs remain fixed. By default, the Yukawas come from the historical A(-2) running-low fitted configuration. An explicit second option injects the frozen-low overnight A(-2) Yukawas as a **hypothetical UV Yukawa input**, not a fitted UV solution.

It is a feasibility test, not a new neutrino fit, and it does **not** claim to rule out UV-running T3 models if no positive UV scalar boundary exists for these fixed Yukawas.

## Numerical calculation

At fixed Yukawas and other initial inputs, the one-loop beta functions for the two scalar mass-squared parameters form an inhomogeneous linear system, while the dimensionless RGEs and fermion mass evolution are independent of these scalar squared masses in the examined truncation. Their combined UV+EFT1 trajectory therefore defines a candidate affine map

\[
\begin{pmatrix} m_{S_1}^2 \\ m_{S_2}^2 \end{pmatrix}_{\mu_S=10^3\,\mathrm{GeV}}
= A\begin{pmatrix}m_{S_1}^2\\m_{S_2}^2\end{pmatrix}_{\Lambda=10^7\,\mathrm{GeV}}+b,
\]

with input/output mass squares expressed in **TeV\(^2\)**. The script does not assume a specific numerical `A,b`: it obtains them from the *real repository's numerical ODE solvers*. It runs the historical UV boundary and independent perturbations of `4 TeV^2` (default), then tests the affine formula against **two extra full RG trajectories** (zero UV masses and a simultaneous two-mass perturbation). The calculation fails closed if the error exceeds `0.05 TeV^2`.

It then:

1. Solves algebraically for UV masses that would give **2 TeV scalars at muS**, retaining negative `m^2` solutions *only as diagnostics*.
2. Uses a constrained linear program to determine whether any **nonnegative UV mass-squared inputs** can place *both* low-scale scalar masses in the fixed 0.5–5 TeV range.
3. Independently checks the more restrictive condition that UV scalar masses themselves lie in 0.5–5 TeV.
4. If feasible, reruns the **full UV+EFT1 integration** at each proposed positive UV input, and checks the original mass hierarchy, scalar threshold and perturbativity screens. An affine feasibility result does not automatically become a physically screened accepted input.

The original fixed matching scales remain:

- `muUV = 1e7 GeV`
- `muF = 1e5 GeV`
- `muS = 1e3 GeV`

The target 0.5–5 TeV band, `rhoF <= 0.5` and Yukawa/quartic perturbativity are **not relaxed**.

## Full files in this package

- `scripts/ProbeT3ReducedUVScalarBoundaries.py` — complete new pilot source
- `tests/python/test_t3_reduced_uv_scalar_boundary_pilot.py` — offline tests for the affine reconstruction, numerical validation, signed roots and constrained mass-window feasibility
- `docs/T3_REDUCED_UV_SCALAR_BOUNDARY_PILOT.md` — this procedure

All outputs are isolated under `output/studies/reduced_scalar_uv_boundary_pilot/<run-id>/`; each run ID is write-once and the script refuses an existing run directory. No previous frozen fit, UV campaign or GitHub source is modified.

## Windows PowerShell

Extract the ZIP into your local checkout:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath ".\T3_Reduced_UV_Scalar_Boundary_Pilot.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_reduced_uv_scalar_boundary_pilot.py
```

Run the baseline A(-2) pilot:

```powershell
python -u -m scripts.ProbeT3ReducedUVScalarBoundaries --run-id a_m2_historical_v1
Get-Content -Encoding UTF8 output\studies\reduced_scalar_uv_boundary_pilot\a_m2_historical_v1\summary.md
```

Optional independent Yukawa-boundary sensitivity check, using the A(-2) best fit from the **completed frozen-low overnight run**:

```powershell
python -u -m scripts.ProbeT3ReducedUVScalarBoundaries --run-id a_m2_frozen_yukawas_v1 --yukawa-source frozen-newfit
Get-Content -Encoding UTF8 output\studies\reduced_scalar_uv_boundary_pilot\a_m2_frozen_yukawas_v1\summary.md
```

Output per run:

- `result.json`: exact forward trajectories, scalar mass² affine map, validation errors, inferred signed target and LP feasibility, physically screened follow-up trials, and input source digests.
- `summary.md`: concise result and caveats.

The script reuses `scripts/RunT3GeneralScalar48Campaign.py`, `scripts/RunT3ReducedScalarSobolComparison.py`, the existing Stage31 UV RGEs and the prior native EFT1 beta export. If the native EFT1 export is not already present, the previous reduced-scalar runner may invoke local `wolframscript` to generate it in this new isolated output directory.

## Interpretation

If both constrained programs are infeasible, that means **no choice of the two specified positive UV squared-mass parameters** achieves two light scalars at the lower threshold **for the chosen fixed Yukawas and other inputs**, in the validated candidate one-loop affine approximation. This is not a proof against simultaneously changing Yukawas, heavy fermion masses, quartics, thresholds or allowing different physics.

If a negative UV mass-squared input would be necessary for a 2 TeV output, it violates the selected positive-UV benchmark assumption; the script will not use that point in a physical running calculation.

If positive UV mass-squared solutions exist, they must also pass the original scalar/fermion hierarchy and perturbativity checks before being considered a viable next-stage UV fitting boundary. No vacuum stability, positive masses *throughout* the whole RG trajectory, or independent physical Matchete/RGBeta coupling mapping is certified.
