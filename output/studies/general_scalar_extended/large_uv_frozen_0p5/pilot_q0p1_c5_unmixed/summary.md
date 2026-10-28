# B/C historical-input general scalar paired pilot

Study: `large_uv_frozen_0p5`; run ID: `pilot_q0p1_c5_unmixed`

Candidate native RGBeta UV/EFT1; distinct zero-extra versus nonzero-extra scalar boundaries.
The nonzero extras are independent benchmark assumptions, not fits or physical coupling normalisations.
Original fit quality is unchanged; no Yukawa refit. No general physical C5 or vacuum-stability result.

| Model | Original fit chi2 | Pass <=5.86 | Shift abs(lambdaT3) (fraction) | Shift mS1² | Shift mS2² | Direct C5 preview (uncertified) |
|---|---:|---|---:|---:|---:|---:|
| B_-1 | 32.2945 | NO | 0.000549799 | 0.000632149 | 0.000833727 | 0.779235 |
| C_-1 | 4.50018e-17 | yes | 0.000549799 | 0.000632149 | 0.000833727 | 0.982239 |

## RGE-generated mixed quadratic coefficient

Zero quadratic *input* does not guarantee a zero quadratic *endpoint*.
| Model | Restricted endpoint abs(mu2)/min(m²) | Quartics endpoint abs(mu2)/min(m²) |
|---|---:|---:|
| B_-1 | 0 | 0.0205559 |
| C_-1 | 0 | 0.0205559 |

## Interpretation

- This is a controlled change of candidate scalar boundary conditions on each **historical saved fitted input**. The output is NOT an updated or refitted fit.
- The C5 preview, if requested, varies only the two direct quartics at the restricted endpoint mass/Yukawa background, with **mu2 explicitly projected to zero for both evaluations**. It assumes an **uncertified identity** RGBeta-to-Matchete map and uses the old one-scale expression.
- One-loop running may induce mu2 even when its boundary value is zero. Such endpoint mixing is tabulated above and in the C5 preview JSON, but is **omitted** from the hypothetical unmixed one-scale C5 sensitivity. This is not the full RG-endpoint C5 change.
- When mu2 is nonzero, mixing can modify C5 at leading one-loop order. No matching is claimed for the actual mixed-mass endpoint or for the optional explicit mixing scenario.
- Historical input/result locations are never written; each run ID is new-only.
