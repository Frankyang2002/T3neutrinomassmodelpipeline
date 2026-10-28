# Six accepted B/C general-scalar candidate benchmarks

Mode: live RG / conditional one-scale C5 replay

| Study | Model | Saved chi2 | rhoF at 100 TeV | Scalars at 1 TeV (TeV) | UV lambdaT3 | Replay |
|---|---|---:|---:|---|---:|---|
| small_uv_running_5e-9 | B_-1 | 0.0028880425 | 0.028489 | 1.839–2.907 | 8.99752e-06 | CandidateReplayAgrees |
| small_uv_running_5e-9 | C_-1 | 3.1330118e-08 | 0.027951 | 2.23–3.217 | 2.13355e-05 | CandidateReplayAgrees |
| small_uv_frozen_5e-9 | B_-1 | 3.2154541e-12 | 0.011429 | 1.016–1.164 | 5e-09 | CandidateReplayAgrees |
| small_uv_frozen_5e-9 | C_-1 | 2.0454128e-17 | 0.011429 | 1.016–1.164 | 5e-09 | CandidateReplayAgrees |
| large_uv_frozen_0p5 | B_-1 | 0.0329504 | 0.011429 | 1.017–1.165 | 0.5 | CandidateReplayAgrees |
| large_uv_frozen_0p5 | C_-1 | 4.5828604e-11 | 0.011429 | 1.017–1.165 | 0.5 | CandidateReplayAgrees |

## Validation

- All six saved chi2 values <= 5.86; predictions, five-variable NuFIT target, Yukawa configurations and source SHA-256 fingerprints checked.
- Scalar eigenmasses are running 2×2 quadratic eigenmasses, not pole masses.
- UV-running cases have the sampled hierarchy ratio from their individual constrained fits; frozen cases are checked against saved fermion masses.
- The UV-running topology couplings were fitted freely and are NOT equal to the historical fixed 5e-9 boundary.
- All six candidates were independently replayed with the existing native candidate RGBeta + one-scale C5 + SM/Weinberg evaluator.
- Diagnostic coordinate screens: max |lambda| <= 4pi (12.566) and max Yukawa singular value <= sqrt(4pi) (3.5449).
- These coordinate screens are not a general-scalar vacuum-stability proof or a physical unitarity calculation.

| Study | Model | max |lambda| | max Yukawa singular | Quartic/Yukawa screens |
|---|---|---:|---:|---|
| small_uv_running_5e-9 | B_-1 | 0.33894 | 0.0265 | PASS |
| small_uv_running_5e-9 | C_-1 | 0.33729 | 0.025442 | PASS |
| small_uv_frozen_5e-9 | B_-1 | 0.27765 | 1.7245 | PASS |
| small_uv_frozen_5e-9 | C_-1 | 0.27765 | 1.4094 | PASS |
| large_uv_frozen_0p5 | B_-1 | 0.50827 | 0.0001468 | PASS |
| large_uv_frozen_0p5 | C_-1 | 0.50827 | 0.00030471 | PASS |

## Scope boundary

- B(-1)/C(-1) only: this is 6/6 available candidate fits, not 48/48 general-scalar fits.
- The Matchete–RGBeta coupling map is conditional; the saved C5 is one-scale, not fully matched sequential EFT.
- Positive sampled scalar quadratic eigenvalues do not prove boundedness of the full potential.
- Grouped matching at 100 TeV / 1 TeV is approximate even when masses satisfy the screen.
- Historical and fitted inputs are never overwritten.
