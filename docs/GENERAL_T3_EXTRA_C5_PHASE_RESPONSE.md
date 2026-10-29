# General T3: extra-only neutrino mass response and phase analysis

**Scope:** Ordinary T3-B and T3-C at alpha=-1, **two independent scalars S1 and S2**. This is not the physical one-scalar shared-scalar scotogenic theory.

## Verified prior to this batch by the user's local execution

- 16/16 ordinary T3 points passed Matchete scalar-quartic validation and extended one-loop matching.
- Formal zero-extra limit reproduced the original one-loop C5 in all 16 cases.
- For B/C(-1) the matched C5 has two new quartic response terms, with y1*y1 and y2*y2 structures.
- Source betas for the two corresponding RGBeta quartics are nonzero at zero extra RGBeta quartics.
- Matchete tensor factors are exactly 1 relative to the explicit canonical SU(2) tensors, for four model/coupling combinations.
- Stage 11 produced extra-only hard C5 matrices for both models, using the toy benchmark. **RGBeta-to-Matchete normalization is NOT yet verified.**

## What Stage 12 actually calculates

For each model, read the `output/general_t3_numerical_pilot/<model>_matchete_hard_c5.json` output from Stage 11. Reconstruct its two complex symmetric C5 contributions C1 and C2 and check:

    C1 + C2 == delta_C5

Use the existing project mass convention:

    delta_m_nu [eV] = -(246.22 GeV)^2/2 * delta_C5 [GeV^-1] * 1e9.

Calculate the singular values of the complex-symmetric delta_m_nu, which coincide with its Takagi singular values, in increasing order. These are NOT physical neutrino mass predictions, because the original C5 and subsequent RGE/mass mixing have not been included.

Report the norm interference:

    ||C1+C2||_F^2 = ||C1||_F^2 + ||C2||_F^2 + 2 Re Tr(C1^dag C2).

Create the zero-extra, either-extra-only, and both-extra benchmark scenarios. Scan the two **physical quartic coupling** phases over 0, 90, 180, 270 degrees (16 combinations). `Bar[Coupling]` in the matched coefficient means rotating the underlying coupling by exp(+i theta) rotates its contribution by exp(-i theta). The scan reuses fixed one-loop kernels and magnitudes, not re-matched Lagrangians.

The Stage 11 user's log gives C5 Frobenius norms of 2.577247421905154e-17 GeV^-1 for B(-1) and 1.4879744928052096e-17 GeV^-1 for C(-1), corresponding to extra-only mass-matrix Frobenius norms about 7.81e-4 and 4.51e-4 eV.

## Files

- `Numerical/experimental/GeneralT3NeutrinoImpact.py`: checked C5/mass matrices, Takagi singular values, interference, phase scenarios.
- `scripts/RunGeneralT3ExtraC5PhaseResponse.py`: CLI for both models, JSON + CSV.
- `scripts/RunGeneralT3ExtraC5PhaseResponse.ps1`: PowerShell wrapper, reuses previous output.
- `tests/python/test_general_t3_extra_c5_phase_response.py`: regression tests, including correct complex-conjugate phases and strict failures.
- `Numerical/experimental/GeneralT3RGSourcePilot.py`: only converts NumPy scalar output fields to plain Python floats; physics unchanged.

## Commands (project root, PowerShell)

    Expand-Archive -LiteralPath .\T3_General_Scalar_Batch12.zip -DestinationPath . -Force
    python -m pytest -q tests\python\test_general_scalar_trial.py tests\python\test_scalar_quadratic_prototype.py tests\python\test_general_scalar_c5_proof.py tests\python\test_general_rgbeta_coverage.py tests\python\test_general_t3_scalar_quartic_response_audit.py tests\python\test_general_t3_matching_r_g_source_pilots.py tests\python\test_general_t3_extra_c5_phase_response.py
    powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3ExtraC5PhaseResponse.ps1

No Matchete/RGBeta reruns are necessary. Results:

    output/general_t3_neutrino_impact/B_m1_extra_mass_impact.json
    output/general_t3_neutrino_impact/C_m1_extra_mass_impact.json
    output/general_t3_neutrino_impact/extra_mass_impact_summary.json
    output/general_t3_neutrino_impact/extra_mass_impact_summary.csv

## Remaining physics — not implemented

1. Gauge-invariant **and convention-safe** identification of the Matchete extra couplings with the RGBeta couplings. Canonical SU(2) agreement on Matchete's side is insufficient: RGBeta `delS2` normalization and the coupling conventions require explicit checks. RGBeta documentation states that `AddQuartic` assumes no normalization factors; the supplied `GroupInvariant` fixes them (Thomsen, EPJC 81 (2021) 408, section 4.2).
2. Complete RGEs for all additional quartics and scalar quadratic mass terms.
3. Extended running states, heavy/scalar mass-basis consistency, and hierarchical matching with all nonzero UV terms.
4. Comparison of the complete low-scale C5 and neutrino masses against the user's fitted benchmark.

**Do not pass Stage 12 delta_C5 as the full physical C5 to the production neutrino solver.**
