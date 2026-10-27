# T3 B_-1/C_-1: direct one-scale C5 to low-scale SM+Weinberg endpoint

## Scope

This is **not** the fermion-first EFT1 matching calculation. It continues the
**already computed** one-scale, all-heavy Matchete C5 for B_-1/C_-1 to the
frozen reference `mu_low_gev` using the project's existing
`Numerical.running.SMWeinbergEvolution.evolve_sm_weinberg` solver.

The method is deliberately a *SM-only-background diagnostic*:

1. Read immutable `Numerical/benchmarks/fitted_reference/benchmark_manifest.json`
   and the model's exact frozen reference JSON, enforcing its SHA256 hash.
2. Read the existing scalar-aligned direct-C5 JSON from
   `output/general_t3_numerical_c5_frozen_scalar_scenario/` and refuse
   anything marked physical or sequential EFT1. Its scalar diagonals come
   from the frozen input, with a separately declared quadratic mixing scenario.
3. Evolve the frozen **SM** renormalisable parameters from its stored UV scale
   (benchmark 10^7 GeV) to the stored scalar threshold (benchmark 10^3 GeV)
   with `K=0`, deliberately omitting heavy/BSM effects in this segment.
4. Insert the **unaltered** one-scale Matchete-derived symmetric C5 at the
   scalar threshold and run the existing one-loop SM+Weinberg equations to the
   saved low scale (benchmark 100 GeV).
5. Calculate the Takagi masses as singular values of `m_nu = -(v^2/2) C5`
   with the frozen benchmark's constant `vev_gev`.

The SM-only evolution above the scalar threshold is **not the actual T3
renormalisable running**. No F threshold, EFT1 LLSS, EFT1 direct Weinberg
running, scalar hard threshold or Matchete/RGBeta nonzero-quartic normalization
has been computed here. It does **not** reproduce the original fitted
oscillation result and must not be labelled as such.

## Project files

- `scripts/RunGeneralT3LowScaleOneScale.py` (new)
- `tests/python/test_general_t3_low_scale_one_scale.py` (new)
- `docs/T3_GENERAL_T3_ONE_SCALE_SM_WEINBERG_ENDPOINT.md` (new)

No existing source or previously generated results are overwritten. The runner
imports production `SMWeinbergEvolution.py` from the user's installed repo;
its high-level regression tests inject a deterministic evolution double.
Those tests are **not** independent verification of the external numerical RGE.

## Windows PowerShell

From repository root:

```powershell
Expand-Archive -LiteralPath "C:\path\to\T3_General_OneScale_SM_Weinberg_Endpoint.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_low_scale_one_scale.py
python -m scripts.RunGeneralT3LowScaleOneScale --all
```

Prerequisite:

```powershell
python -m scripts.RunGeneralT3FrozenScalarC5 --all
```

Expected outputs:

- `output/general_t3_low_scale_one_scale/B_m1_final_low_scale_one_scale.json`
- `output/general_t3_low_scale_one_scale/C_m1_final_low_scale_one_scale.json`
- `output/general_t3_low_scale_one_scale/summary_final.json`

The JSON reports scalar-boundary and low-scale Takagi singular masses, change
ratios, C5 matrices, source/reference hashes, and explicit limitations.

## Remaining physical step

The proper fermion-first T3 calculation requires the mass-eigenfield
Lagrangian with F integrated out at `mu_F` while X,Y remain active,
the intermediate LLSS/direct Weinberg running with the appropriate tensor
basis, and the scalar matching at `mu_S` with correctly mapped couplings.
Those results must replace—not be conflated with—this diagnostic.
