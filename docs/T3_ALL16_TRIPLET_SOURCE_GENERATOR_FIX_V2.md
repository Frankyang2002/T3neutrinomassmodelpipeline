# All-16 reduced-potential full Wilson trial v2 — native triplet source basis

## Scope

A pure basis correction to the existing all-16 exploratory full-Wilson trial.
The previous `6/16` run is preserved; its blocked A/D/E results are not physical failures.
No production RGBeta source, fermion threshold, scalar matcher or GitHub content is modified.
The existing `scripts/RunT3All16ConditionalFullWilsonTrial.py` CLI is unchanged.

## Verified root cause

A/D/E contain SU(2) triplet scalar representations whose *native Matchete
Cartesian components* are not the spin-1 weight (`m=+1,0,-1`) components
used by `Aminus4EFT1Couplings.su2_generators(3)`.

With the Matchete source-field triplet `S_native` and spherical field `S_weight`:

```
S_weight = N @ S_native
N = [[-1/sqrt(2), 0, 1j/sqrt(2)],
     [0,           1, 0          ],
     [ 1/sqrt(2), 0, 1j/sqrt(2)]]
T_native^A = N.conj().T @ T_weight^A @ N
```

This is the same SU(2) intertwiner determined from the existing A(-4) Matchete
matched CG source. It transforms Hermitian spin-one generators into the
pure-imaginary antisymmetric Cartesian matrices required by the release.
The original v1 trial instead used `T_weight^A` on native triplet components.

Independent algebra checks using *actual* nonzero `p=q=1` Matchete LLSS
source components from A(-4), D(-2), E(-4):

| Model | Original wrong spherical generator: max SU(2) Ward residual | Correct Cartesian triplet generator |
|---|---:|---:|
| A_-4 | 2.0 | < 3e-16 |
| D_-2 | 2.0 | < 3e-16 |
| E_-4 | 2.0 | < 3e-16 |

This explains the ~0.023–0.026 full Weinberg tensor projection residual
after running for the ten A/D/E models, *without assuming they have zero C5*.
It has not yet been demonstrated that all ten previously blocked integrated
runs succeed in the local Windows checkout.

## Changes in ZIP (full replacements)

- `Numerical/running/All16ConditionalFullWilsonTrial.py` — full replacement.
  1. Use the native Cartesian triplet SU(2) generator on any scalar triplet,
     and transform generators with the existing mixed/split real rotation.
  2. Before each Wilson ODE, calculate source and reduced-quartic SU(2)/U(1)
     Ward residuals in the *actual matching basis*; stop on a real violation.
  3. Include Ward residuals in per-model JSON and comparison CSV.
  4. Keep the six-term Wilson beta, source normalization convention, RK2
     integrator, reduced potential approximation, and joint scalar matcher.
- `tests/test_all16_conditional_full_wilson_trial.py` — full replacement with
  12 tests, including three source-component witnesses from Matchete.
- `docs/T3_ALL16_TRIPLET_SOURCE_GENERATOR_FIX_V2.md` — this file.

## PowerShell

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive -LiteralPath "$HOME\Downloads\T3_All16_Triplet_Cartesian_Ward_Fix_v2.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p 'test_all16_conditional_full_wilson_trial.py' -v

# Use fresh destinations: the v1 runs already occupy their directories.
python scripts\RunT3All16ConditionalFullWilsonTrial.py --model A_-4 --steps 2 `
  --output-root 'output\general_observable\all16_full_wilson_A_m4_smoke_v2'
python scripts\RunT3All16ConditionalFullWilsonTrial.py --steps 2 `
  --output-root 'output\general_observable\all16_full_wilson_reduced_potential_v2'
```

Inspect `output/general_observable/all16_full_wilson_reduced_potential_v2/summary.json`
and `all16_trial_comparison.csv`. Each model has `source_Ward_max`,
`quartic_Ward_max`, `max_weinberg_residual`, mass estimates, or an explicit
blocker. Failures are not silently projected away.

## Numerical approximation retained

This remains a **rough conditional** full-Wilson trial: the earlier v6 scalar
background supplies only `HHSS`, and the transport includes `lambdaH` plus
that `HHSS` tensor (other quartics omitted). It linearly interpolates SM
couplings between the saved F/S scales and freezes the HHSS tensor at S.
The 16-model calculation is therefore **not** equivalent to A(-4)'s full
17-beta corrected RGBeta + nine-quartic trajectory. Neither Matchete signed
CC nor hard fermion one-loop matching is inferred.

EFT order is unchanged:
`UV -> integrate out ALL heavy fermions -> EFT1 (S1 and S2) -> integrate out ALL heavy scalars simultaneously -> SMEFT`.
