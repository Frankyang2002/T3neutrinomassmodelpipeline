# A(-4) 18-real-Yukawa rough fit v1

This is a deliberately approximate numerical continuation of the successfully
executed A(-4) direct hard/soft finite channel and SMEFT/PMNS benchmark. It
fits the **18 real UV y1/y2 entries** to the five NuFIT 6.1 NO observables.

## Installed project files

- `Numerical/parallel_t3/complete/a4_rough_18fit.py` — full numerical fit implementation (new).
- `scripts/RunT3A4Rough18RealFit.py` — command-line entry point (new).
- `tests/test_t3_a4_rough_18fit.py` — unit tests (new).
- `docs/T3_A4_ROUGH_18_REAL_FIT_v1.md` — this note (new).

The existing first-threshold and scalar-matching WXF/Matchete calculations are
**not** repeated or modified. The previous `a4_direct_hard_soft.py` and
`a4_rough_benchmark.py` modules must already be installed.

## Fixed physics assumptions and calculation

- UV 1e7 GeV -> F 1e5 GeV -> joint (S1,S2) 1e3 GeV -> SMEFT 100 GeV.
- Both Yukawas are 3x3 **real at UV**; all 18 entries are optimizer variables.
- The baseline UV-to-F response of each Yukawa entry is fixed to
  `yF_baseline[p,r] / yUV_baseline[p,r]`. The background masses, quartics, and
  SM RG coefficients are held at the saved benchmark. This does NOT perform
  complete UV RG running at each trial.
- At muF, the source-grounded ordered direct coefficient is computed via
  `formal_numeric_direct_remainder`; the provisional symmetric coefficient is
  `K(muS)=A+A.T`, with GammaCC/Weyl factor assumed 1.
- No separately computed EFT1 tree-to-LLHH logarithm is added, avoiding
  double counting of the shared direct hard/soft logarithm.
- The existing `SMWeinbergEvolution.evolve_sm_weinberg` provides six basis
  evolutions of a symmetric complex 3x3 `K`; this yields an exact complex-linear
  map `K(muS) -> K(100 GeV)` with the saved SM running background. All optimizer
  evaluations use that precomputed map; the final point is recomputed through
  the existing SMEFT solver and compared to the map.
- The existing `rough_numerics.neutrino_observables` supplies Takagi/PMNS; the
  five Gaussian residuals use `configs/nufit_6_1_ic24_no.json`.
- Initial guess scales both Yukawa matrices by
  `(NuFIT_delta_m31_sq / baseline_delta_m31_sq)**(1/4)`, approximately 0.74705.
- `scipy.optimize.least_squares` fits 18 real parameters; a small regularization
  about that starting guess removes gross parameter drift in this underdetermined
  five-observable/18-parameter system. Reported chi2 excludes regularization.

## Local commands (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_A4_Rough_18Real_Fit_v1.zip" -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p test_t3_a4_rough_18fit.py -v
python -X utf8 -m scripts.RunT3A4Rough18RealFit
Get-Content "output\parallel_t3\a4_rough_18real_fit_v1\A_-4\fit.json"
```

`fit.json` provides starting/fitted chi2, low-energy masses, mixing angles,
18 fitted UV values in units of `1e-4`, and the approximate F Yukawas.
`fitted_matrices.npz` stores the corresponding 3x3 complex matrices.

For a different optimizer budget without overwriting results:

```powershell
python -X utf8 -m scripts.RunT3A4Rough18RealFit --max-nfev 400 --out output\parallel_t3\a4_rough_18real_fit_v1_more\A_-4
```

## What this result does NOT mean

This is an approximate fit to the direct hard/scalar-soft channel, **not** a
complete finite matched physical `C5`. It does not recompute full UV Yukawa RG
at each point, include every first-threshold dimension-five Wilson operator,
run all EFT1 Wilson operators consistently with the counterterm subtraction,
or include the remaining finite two-scalar matching terms. The 18 real fitted
UV entries are parameters of the stated *fixed-background transfer model*.

The objective is to get useful rough masses/angles and a fitted benchmark now,
without expanding the project into more certification or diagnostics.
