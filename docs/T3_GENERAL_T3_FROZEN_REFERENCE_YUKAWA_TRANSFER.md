# General T3 one-scale frozen-reference Yukawa transfer (B_-1/C_-1)

## Purpose and limits

This is **NOT** a physical fitted-reference re-evaluation, a new fit, or
sequential fermion-first EFT1 matching. It is a controlled *hybrid* comparison:

- **Yukawas and fermion masses** are taken from the existing frozen fitted
  reference files in `Numerical/benchmarks/fitted_reference/`;
- **scalar masses and mixing U** remain those from the independently validated
  `output/general_t3_mass_eigen_matchete/{B,C}_m1_final_parameters.json`;
- **one-loop amplitude** is evaluated from the original Matchete InputForm
  coefficient, not a phenomenological substitute loop function;
- **normalisation** uses the project's ordered `A5` -> symmetric `C5`
  convention; full Matchete/RGBeta component map is *not* certified.

The frozen reference scalar parameters differ from the experimental scalar
mass-mixing inputs, even before considering the extra off-diagonal mass term.
Consequently, **the resulting neutrino masses are not the original fitted
NuFIT predictions**, and cannot be used to claim a successful fit or
reproduce the legacy `smallY_smallL` benchmark. The output includes the
original and experimental squared scalar mass values and their differences.

The two additional quartic couplings in the exported coefficients are set to
zero only because their corresponding frozen-reference input couplings are
**explicitly zero**. If either input is nonzero, the runner fails closed:
its Matchete/RGBeta operator-normalisation map is not yet certified.

## Flavour and mass treatment

The frozen `ordinary.y1`, `ordinary.y2` are 3x3 complex matrices with three
fermion-flavour columns. `ordinary.MF` is a *diagonal*, real-positive 3x3
mass matrix. The evaluator computes the scalar expression independently for
each fermion-flavour column `k`, substituting `MF = MF[k,k]` in that term:

`A5[i,j] = sum_k A5_ Matchete [y1[i,k], y2[j,k], MF[k,k], MSX, MSY, U, quartics]`.

This implementation assumes no fermion mass rotation, no off-diagonal heavy
fermion mass, and the same scalar eigenstate benchmark for all three fermion
flavours. Off-diagonal/complex `MF` fails closed rather than being treated as
a scalar number. The original one-column illustrative mode still works.

The `mnu` matrix uses `C5 = A5 + Transpose[A5]` and
`mnu = -v^2 C5 / 2`, with `hbar = 1/(16 pi^2)` once. The three Takagi masses
are the singular values of the complex-symmetric mass matrix, in eV.

`max_abs_mnu_eV` is **the largest absolute matrix entry** `max_ij |mnu_ij|`.
It is *not* the largest physical neutrino mass. For the previous illustrative
inputs, the largest Takagi mass is 148.42297 eV (B_-1) or 86.02098 eV
(C_-1); the largest matrix entries are only 61.20390 and 35.47340 eV.
The updated CLI reports both quantities by their unambiguous names.

## Installation

Extract the ZIP at the repository root. Changed file paths:

- `scripts/RunGeneralT3NumericalC5.py` — FULL replacement; supports a
  nondegenerate `fermion_masses_gev` list as well as the original common MF.
- `tests/python/test_general_t3_numerical_c5.py` — FULL replacement (retains
  existing tests).
- `scripts/RunGeneralT3ReferenceC5.py` — NEW frozen-reference transfer CLI.
- `tests/python/test_general_t3_reference_c5.py` — NEW regression tests.
- `docs/T3_GENERAL_T3_FROZEN_REFERENCE_YUKAWA_TRANSFER.md` — NEW note.

Frozen reference configurations are **read-only**. Nothing is written to their
original files, RGBeta exports, or Matchete model builder.

### Commands (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -m pytest -q tests\python\test_general_t3_reference_c5.py tests\python\test_general_t3_numerical_c5.py
python -m scripts.RunGeneralT3ReferenceC5 --all
```

To reproduce the earlier illustrative run with clarified labels:

```powershell
python -m scripts.RunGeneralT3NumericalC5 --all
```

The transfer reads:

- `Numerical/benchmarks/fitted_reference/benchmark_manifest.json`
- `Numerical/benchmarks/fitted_reference/T3_dS1_2_dS2_2_dF_{1,3}_alpha_m1.json`
- `output/general_t3_mass_eigen_matchete/{B,C}_m1_final_parameters.json`
- `output/general_t3_mass_eigen_matchete/one_loop_probe/{B,C}_m1_final_one_loop_probe.json`
- corresponding `*_weinberg_ordered_coefficient.inputform.txt` files.

Outputs (separate from the illustrative set):

`output/general_t3_numerical_c5_reference_transfer/{B,C}_m1_final_numerical_c5.json`

`output/general_t3_numerical_c5_reference_transfer/{B,C}_m1_final_mnu_eV.csv`

`output/general_t3_numerical_c5_reference_transfer/summary_final.json`

Each model's JSON records the frozen input SHA256 and mass mismatch, with
`physical_c5_certified=false`, `sequential_EFT1_performed=false`, and
`mass_basis_compatible_with_frozen_reference=false`.

## Next physical milestone

Do not conflate this one-scale hybrid transfer with the source-consistent
fermion-first EFT1 result. Physical next steps are (i) choose/derive the
quadratic mixing input and mass basis **at the same running scale** as the
EFT1 benchmark, (ii) match the Matchete and RGBeta component-coupling basis,
(iii) run the appropriate F -> scalar thresholds, and (iv) construct the
3x3 symmetric C5 and masses at the specified comparison scale.
