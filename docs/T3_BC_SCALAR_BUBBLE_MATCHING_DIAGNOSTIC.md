# B(−1)/C(−1): finite scalar-bubble matching diagnostic

## Scope and provenance

This experimental package **does not claim a complete physical one-loop scalar-threshold match**. It computes the renormalized, zero-momentum **single insertion of the project's scalar-pair Wilson operator and scalar quartic** at the running Y threshold. It then propagates that *candidate* shift through the existing EFT2(H,X) solver and the SM+Weinberg solver, without changing previous results.

The implementation is derived from the project source conventions:

* `RGE/general/WilsonTensorRGE.py`: the one-loop quartic contribution to the Wilson beta is `+ lambda_abcd C_ijcd` (one-loop numerator).
* `Numerical/experimental/T3BCCoupledWilsonTensor.py`: fully symmetric 12-real-scalar `C_ijabf`, full 4-index quartic potential tensor `lambda_abcd`, and six physical lepton flavour pairs.
* `Numerical/experimental/T3BCScalarYThresholdProjection.py`: `(H,X,Y)` mass-basis rotation, including epsilon and general complex U.
* `Numerical/experimental/T3BCEFT2HX.py`: `H,X` Wilson running and scalar RGEs.
* `scripts/RunT3BCEFT2Sequential.py`: saved Y/X scalar thresholds and conditional low-energy conversion.

No independent Matchete amplitude check of the finite coefficient has been performed. The sign is fixed **relative to the recorded RGE convention** by renormalization-scale invariance, but this alone does not certify the full matching normalization.

## Scalar-pair beta and finite loop

For `V_4=lambda_abcd phi_a phi_b phi_c phi_d/4!`, the project RGE quartic term is

\[
16\pi^2\frac{d C_{ijab,f}}{d\ln\mu}\bigg|_{\lambda}=
\sum_{cd}\lambda_{abcd}C_{ijcd,f}.
\]

The finite zero-momentum scalar two-propagator loop in \(\overline{\mathrm{MS}}\) is

\[
\overline B_0(0;a,b;\mu^2)=
1-\frac{a\ln(a/\mu^2)-b\ln(b/\mu^2)}{a-b}.
\]

It tends to \(-\ln(a/\mu^2)\) when \(b\to a\), and to \(1-\ln(a/\mu^2)\) when \(b\to0\). At fixed masses,

\[
\frac{d\overline B_0}{d\ln\mu}=2.
\]

The **candidate scalar-pair contribution** inferred using this RGE normalization is

\[
\Delta C_{ijab,f}^{\mathrm{pair}}
=-\frac{1}{2(16\pi^2)}\sum_{cd\;:\;c\text{ or }d\in Y}
\lambda_{abcd} C_{ijcd,f}\,\overline B_0(0;m_c^2,m_d^2;\mu_Y^2),
\quad a,b\in\{H,X\}.
\]

A direct differentiation reproduces minus the heavy-leg contribution to the full-theory scalar-pair beta, as required when passing into EFT2. The code checks the full partition into `HY`, `XY`, and `YY` internal blocks and preserves all six flavour pairs.

## Mass hierarchy and two regimes

From the successful user-local conditional EFT2 runs (the mass at the X threshold is used as a scale proxy only),

| | B(−1) | C(−1) |
|---|---:|---:|
| \(\mu_Y\) (GeV) | 1187.805782 | 1185.393345 |
| \(\mu_X\) (GeV) | 1035.777755 | 1032.855235 |
| \(\mu_X/\mu_Y\) | 0.872009 | 0.871319 |

The **exact** Y-boundary mass ratio in the calculation is extracted from the saved Y crossing's *same-scale* \(m_X^2(\mu_Y)\), not from the lower EFT2 crossing \(\mu_X\). These quantities must not be conflated. For ratios of the size quoted above, the representative mixed scalar loop has \(\overline B_0\sim 0.13\) at \(\mu=\mu_Y\), whereas its massless-X limit is precisely `1`.

Both versions are provided:

1. `full-mass`: evaluate the finite zero-momentum scalar bubble using the exact \(m_X^2(\mu_Y)\) and \(m_Y^2(\mu_Y)\). This keeps mass-ratio effects that correspond to higher-dimensional power corrections in a strict EFT expansion. **It is not a certified finite matched EFT.**
2. `hard-massless-X`: set \(m_X^2=0\) *inside the mixed loop*. This gives the formal leading hard-mass-expansion coefficient, but the ratio \(m_X^2/m_Y^2\sim0.76\) is not parametrically small.

The `YY` equal-mass bubble vanishes at \(\mu_Y=m_Y\) **within this one-loop scalar-pair term**, and the corresponding pure `XX` bubble vanishes at \(\mu_X=m_X\). Neither fact implies that the complete finite scalar-threshold match is zero: other topologies, renormalizable threshold counterterms, operators and mixed heavy/light effects remain to be checked.

## Numerical workflow

The runner:

1. Verifies the hash-protected Y projection source and the saved EFT2 baseline.
2. Reads the full `(3,3,12,12,6)` Wilson tensor **already rotated to mass basis** at the Y crossing.
3. Replays the saved full EFT1 RG background and reconstructs the real, fully symmetric quartic tensor at the same crossing. Rotates all four scalar indices using the existing certified field transformation.
4. Checks that the projected 8-real-field quartic tensor reproduces the saved EFT2 six-coupling inert-doublet fit.
5. Constructs both scalar-pair finite candidate tensors, retaining `(3,3,8,8,6)` and all lepton flavour pairs.
6. Runs **both corrected candidates** and an **unmodified control** through the same EFT2 Wilson and SM+Weinberg RG solvers. The unmodified control must reproduce the saved baseline mass singular values to a strict tolerance.
7. Writes separate JSON/NPZ outputs and summary CSV, with explicit physical-certification flags **false**.

## Setup / exact commands (PowerShell)

Install the ZIP at the repository root while keeping the paths within it:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\t3_bc_scalar_bubble_matching_diagnostic_v1.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p "test_t3_bc_scalar_bubble_threshold.py" -v
python -u -X utf8 -m scripts.RunT3BCScalarBubbleCandidate `
  --all `
  --run-id continuous_v1 `
  --y-output-id bc_y_threshold_projection_v1 `
  --base-output-id bc_eft2_sequential_v1 `
  --output-id bc_scalar_bubble_candidate_v1
```

Saved at:

`output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/bc_scalar_bubble_candidate_v1/`

The report includes `candidate_masses_ev_descending`, correction tensor norms, `B0_YX`, change in \(C_5(\mu_X)\), and strict no-correction replay parity. **Neither candidate is a physically certified neutrino spectrum.**

## Remaining tasks for physical completion

* Derive independent **native Matchete** one-loop matching with Y heavy/X light, including the correct subtraction of EFT2 loops and finite counterterms. The available two-scalar Matchete probe integrates X and Y *together*; it is not already that sequential calculation.
* Validate the scalar-pair finite threshold term's **overall normalization and sign** against the native amplitude/Weinberg projector; validate light-heavy mass dependence consistently with EFT2 operator power counting.
* Derive the corresponding X matching and renormalizable quadratic/quartic finite thresholds, including higher-dimensional effects if required by the close mass scales.
* Replace conditional extra-quartic maps and one-generation SM Yukawa approximation with certified full-flavour results before claiming complete phenomenology.
* Preserve `docs/T3_GENERAL_SCALAR_MATCHING_RESOLUTION_LOG.md`: append the included milestone; never change existing physical certification flags.

External context: Liao and Ma, *One-loop matching of scotogenic model onto SMEFT up to dimension 7*, JHEP 12 (2022) 053, arXiv:2210.04270, demonstrates why higher-dimensional effects and explicit one-loop matching matter. This reference is **not** a certificate for the model-specific formula implemented here.
