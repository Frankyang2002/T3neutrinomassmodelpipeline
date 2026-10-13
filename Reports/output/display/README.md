# T3 numerical study — comparison display

This is a compact presentation copy of the standard comparison-only T3 study.
No numerical calculation has been rerun and the full output is unchanged.
Legacy optimal/optimizer outputs are not read or copied.

## Included

- Interactive cross-model comparison dashboards.
- Selected figures for each model, including:
  - UV Yukawa / lambda_T3 running
  - direct intermediate LLSS -> Weinberg generation
  - hard/direct/combined C5 at the scalar threshold
  - final Weinberg coefficient C5 running
  - neutrino mass-splitting and mixing-angle running
- `model_summary.csv` with benchmark metadata and low-energy observables.

## Operator-running convention

The intermediate curve is the implemented direct O(hbar) LLSS -> Weinberg
contribution Delta C5^direct(mu). The final curve is C5(mu) in the
SM+Weinberg EFT. A numerical LLSS self-running trajectory is not inferred
or fabricated by the plotting layer. Feeding one-loop LLSS self-running
through the scalar loop would first modify C5 at O(hbar^2).

## Common comparison benchmarks

The listed lambda_T3 values are common reference values, not necessarily
the coupling initialized in every representation. The comparison uses
lambda_T3^eff = f_R lambda_T3^ref with
f_R(A,B,C,D,E) = (1, 1, sqrt(3), 1, 1/sqrt(2)).

- `smallY_smallL`: Y = 0.005, lambda_T3^ref = 0.1
- `smallY_largeL`: Y = 0.005, lambda_T3^ref = 1
- `largeY_smallL`: Y = 0.5, lambda_T3^ref = 0.1
- `largeY_largeL`: Y = 0.5, lambda_T3^ref = 1

The physical matching/RGE group factors are unchanged; this normalization
acts only on benchmark initialization. `model_summary.csv` records both
lambda_T3^ref and lambda_T3^eff for every successful row.

The comparison mode uses the same fixed non-diagonal Yukawa textures for
every model, scaled by the common Y value. These are controlled comparison
benchmarks rather than oscillation-data fits.

Report model directories copied: 64.
Successful numerical rows summarized: 64.
