# General-scalar B/C pilot against the three historical studies

## Scope

A **separate** exploratory comparison using original saved fitted B(-1)/C(-1)
Yukawa, fermion-mass and scalar inputs. The program reads each source under:

- `output/studies/lambda_5e-9_refit/{fitted,fits}/` — small \(\lambda_{T3}=5\times10^{-9}\), **UV RG on**;
- `output/studies/fermion_start_lambda_comparison/5e-9/{fitted,fits}/` — same topology quartic, **UV frozen**, boundary at fermion threshold;
- `output/studies/fermion_start_lambda_comparison/0p5/{fitted,fits}/` — \(\lambda_{T3}=0.5\), **UV frozen**, boundary at fermion threshold.

No inputs, fit files, historical output, Stage 36, or Batch 7–8 report is modified.
The generated study output is newly created only under
`output/studies/general_scalar_extended/<case>/<run_id>/`. An existing run ID
is rejected; use a distinct `--run-id` for any new scenario or repeated run.

## Benchmarks and actual calculations

For each of B(-1) and C(-1) the operator registry has **seven** additional
complex quartic slots (six existing RGBeta coordinates and one generated
sidecar) and one allowed mixed quadratic operator.

Use the original fitted ordinary scalar boundary, gauge/Yukawa/fermion values
unchanged. At the declared input boundary:

| Scenario | 7 additional quartics | mixed quadratic coefficient |
|---|---|---|
| `restricted_zero` | all exactly 0 | 0 |
| `general_quartics` | each `+0.1 abs(lambdaT3)` by default (real) | 0 |
| `general_quartics_and_mixing` (`--include-mixing`) | same nonzero quartics | `+0.01 min(mS1Sq,mS2Sq)` by default (real, GeV²) |

These are **chosen illustrative candidate RGBeta-coordinate initial values**,
not fitted values or established physical normalizations. These choices do not
impose positivity/vacuum stability on the full quartic potential.

For UV-on: integrate the **full coupled** candidate RGBeta system from 10^7 to
10^5 GeV, tree-project through fermion decoupling, and integrate the full
candidate EFT1 system from 10^5 to 10^3 GeV. For UV-frozen: do **not** run the
UV segment; impose the original fitted input values plus the chosen additional
couplings at 10^5 GeV, project fermions, and run only EFT1 to 10^3 GeV. No
fermion-threshold finite matching of scalar couplings is included.

Compare diagonal running scalar mass parameters, the topology quartic,
all seven additional quartics and (optionally) the scalar-mass eigenvalues and
mixing parameter from the specified SU(2) candidate mass matrix.

## One-scale C5 preview (optional)

`--c5-preview` reads the existing B/C one-scale Matchete ordered LLHH
coefficient and evaluates the symmetrized 3x3 C5 twice. **All** of the
Yukawas, fermion masses, scalar masses, scalar mixing matrix (mu2=0),
\(\lambda_{T3}\) and Matchete symbolic source are held fixed at the
*restricted* scenario values **in an explicitly hypothetical unmixed projection**. Only the two direct additional quartics
are varied to the `general_quartics` EFT1 endpoint values. The diagnostic is

\[
  \epsilon_{C_5}^{\rm identity\ preview} =
   \frac{\|C_5(\lambda_{a,b}^{\rm RGBeta})-C_5(0,0)\|_F}
        {\|C_5(0,0)\|_F}.
\]

**Crucial distinction:** even if the added quadratic coupling is **zero at the
input boundary** of `general_quartics`, coupled one-loop RGEs can induce a
**nonzero quadratic mixing at the EFT1 endpoint**. This is not a failed RG
run. The `--c5-preview` calculation explicitly sets mu2=0 in *both* C5
evaluations to isolate direct quartic insertion at a common unmixed background.
The actual endpoint quadratic coefficients and |mu2|/min(m²) ratios for both
scenarios are reported in the per-model JSON and the summary. No physical C5
for the mixed-mass endpoint is computed; therefore the preview must not be
interpreted as the total C5 difference including mixing.

**Other limitations:** this assumes the **unverified** identity mapping of the
nonzero RGBeta quartics to the Matchete symbols, freezes original study
Yukawas (not UV-evolved Yukawas), and substitutes a simultaneous one-scale
coefficient rather than performing sequential EFT1 Wilson matching. This is
**not** a physical \(C_5\) error estimate, neutrino mass prediction, or
complete general-scalar matching. The *actual* nonzero quadratic mixing is not inserted into the one-scale
\(C_5\) evaluation: it can affect leading one-loop matching and needs a
certified mass-basis/field-map treatment.

The historical fit chi² and `Success`/`PoorFit` flags are only provenance;
new scalar scenarios are not refitted. Do not label a poor historical fit as
physically accepted because the numerical RG solver succeeds.

## Installation and PowerShell execution

Extract the ZIP at the **project root**, keeping the three historical study
folders unchanged. The archive contains **full new files**, no code patches.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

# Test without rerunning the fitted studies:
python -m pytest -q tests\python\test_t3_general_scalar_bc_pilot.py

# Fast pilot: B(-1) with the small UV-frozen study only:
python -m scripts.RunT3GeneralScalarBCPilot --case small_uv_frozen_5e-9 --model B_-1 --run-id pilot_q0p1

# All three historical cases, both models; baseline versus 7 extra quartics:
python -m scripts.RunT3GeneralScalarBCPilot --all --c5-preview --run-id pilot_q0p1_c5_unmixed_projection

# OPTIONAL: mixed-mass scenario also, as a separate NEW result ID:
python -m scripts.RunT3GeneralScalarBCPilot --all --include-mixing --run-id pilot_q0p1_mu0p01

# Read the six summaries:
Get-ChildItem .\output\studies\general_scalar_extended\*\pilot_q0p1_c5\summary.md | ForEach-Object {
    Get-Content -Encoding UTF8 $_.FullName
}
```

Default extra-quartic fraction is 0.1 of the original fitted |lambdaT3|;
`--quartic-fraction 0.01` selects 1%. `--mu2-fraction 0.01` sets 1% of the
smaller initial diagonal mass squared. Both use **positive real** coefficients.
For distinct choices use distinct `--run-id` values.

## Verification status

Offline tests exercise input/scenario validation, path exclusivity, candidate
mass spectra, RG endpoint comparisons and C5 preview guards. Full numerical
integration on the user's Windows checkout is required before quoting any
actual pilot values. The implementation reuses the saved candidate RGBeta
beta exports; no Wolfram regeneration or GitHub changes are necessary.

## Windows zero-input-mixing preview correction

The first runner incorrectly required `quadratic_mu2_gev2 == 0` at the **EFT1
endpoint** to run the optional `--c5-preview`. A vanishing boundary coefficient
can become nonzero under the coupled RGEs. The corrected runner does not erase
or falsify the generated coefficient: it records it in each RG scenario and in
`C5_preview.unmixed_projection_info`, while the *separate direct-quartic*
preview explicitly evaluates a hypothetical diagonal scalar mass matrix.
This does not provide a new full-mixing C5 calculation.
