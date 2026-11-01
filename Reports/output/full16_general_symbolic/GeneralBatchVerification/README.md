# T3: full general-interaction 16-model physics and group-factor batch

**Structural outputs checked:** 16/16.
**Exact RGBeta gauge-group coefficients verified:** 1196/1196. Mismatched: 0, unprojectable: 0.

This is a **read-only** comparison of the already saved and source-fingerprinted general BSM Lagrangians and native RGBeta UV/EFT1 beta functions.

## Deliverables

- `group_factor_coefficients.csv`: exact gauge and universal gauge-linear factors extracted from native beta functions, with independent predictions and residuals.
- `scalar_gauge_structure.csv`: SU(2) and U(1) group-factor predictions for each tested scalar coupling.
- `physical_coupling_dictionary.csv`: 50 extra quartic coupling slots, eight extra bilinears and 38 Yukawa slots, with the actual RGBeta basis block(s) and uncertified physical maps.
- `native_beta_term_index.csv`: **all** native UV/EFT1 beta expressions split into additive terms for inspection. Non-gauge coefficients are indexed, **not independently certified**.
- `model_status.csv`, `summary.json`: machine-readable status and source-qualified audit outcomes.
- `models/`: consolidated per-model overview and links to complete native BSM Lagrangians and all UV/EFT1 RGEs.

## Verified versus outstanding

**Verified by this batch where status says VerifiedRGBetaGaugeCoefficient:** the native RGBeta one-loop coefficient of the indicated gauge monomial agrees exactly with its independent representation formula.

**NOT certified:** the Matchete-to-RGBeta physical action dictionary, every non-gauge Yukawa/quartic group-factor coefficient, dimension-five Wilson RGEs, mass-eigenbasis threshold matching, and C5.

The legacy `Reports/GroupFactors` output uses analytic formulas for the restricted operator basis and cannot automatically certify all newly allowed scalar invariants.

## Per-model status

| Model | Structural | Gauge checks | Matched | Unprojected | Details |
|---|---|---:|---:|---:|---|
| A_+0 | VerifiedStructural | 80 | 80 | 0 | [Report](models/A_p0.md) |
| A_-2 | VerifiedStructural | 80 | 80 | 0 | [Report](models/A_m2.md) |
| A_-4 | VerifiedStructural | 54 | 54 | 0 | [Report](models/A_m4.md) |
| B_+1 | VerifiedStructural | 62 | 62 | 0 | [Report](models/B_p1.md) |
| B_-1 | VerifiedStructural | 90 | 90 | 0 | [Report](models/B_m1.md) |
| B_-3 | VerifiedStructural | 62 | 62 | 0 | [Report](models/B_m3.md) |
| C_+1 | VerifiedStructural | 62 | 62 | 0 | [Report](models/C_p1.md) |
| C_-1 | VerifiedStructural | 90 | 90 | 0 | [Report](models/C_m1.md) |
| C_-3 | VerifiedStructural | 62 | 62 | 0 | [Report](models/C_m3.md) |
| D_+0 | VerifiedStructural | 80 | 80 | 0 | [Report](models/D_p0.md) |
| D_+2 | VerifiedStructural | 54 | 54 | 0 | [Report](models/D_p2.md) |
| D_-2 | VerifiedStructural | 80 | 80 | 0 | [Report](models/D_m2.md) |
| E_+0 | VerifiedStructural | 100 | 100 | 0 | [Report](models/E_p0.md) |
| E_+2 | VerifiedStructural | 70 | 70 | 0 | [Report](models/E_p2.md) |
| E_-2 | VerifiedStructural | 100 | 100 | 0 | [Report](models/E_m2.md) |
| E_-4 | VerifiedStructural | 70 | 70 | 0 | [Report](models/E_m4.md) |

## Next physics checks

1. Reconstruct the physical scalar action for each independent RGBeta quartic tensor, including factorials, H.c., signs, and field-basis conventions, and compare to Matchete component amplitudes.
2. Derive the full signed/complex Yukawa action dictionary from the saved 38 Yukawa CG witnesses, without identifying a CG norm ratio with a physical coupling map.
3. After a certified dictionary exists, compare *all* beta monomial coefficients and repair genuine discrepancies (if found).
