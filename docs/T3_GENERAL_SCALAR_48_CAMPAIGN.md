# 16 × 3 general-scalar campaign — result-recovery fix

## Scope

This is the full replacement for the original 48-study bulk campaign runner,
with **two verified historical-source bookkeeping corrections**. All 16
models, the three scenario labels, the numerical RG equations, the C5 kernels,
and the fitting algorithm remain unchanged.

1. The four B(-1)/C(-1) frozen-UV accepted fits live in files such as
   `B_m1_exploratory_fit.json`, not a `fit.json` at the top of their run folder.
   The fix uses the previously replay-validated six-fit registry at
   `output/studies/general_scalar_extended/validation/bc_six_validation_v1/report.json`
   and verifies the **source fit and fitted-config SHA-256 digests**. All six B/C
   accepted fits are consistently represented in the bulk table, including the
   two full-UV B/C cases whose UV lambdaT3 boundary was relaxed by earlier
   fitting.
2. Original historical **E-class** benchmark source files store the topology
   coupling at one half of the non-E-class value:
   - study `small_uv_frozen_5e-9`: E uses `2.5e-9`
   - study `large_uv_frozen_0p5`: E uses `0.25`
   - study `small_uv_running_5e-9`: E starts with `2.5e-9`
   These are the original stored E benchmark boundaries; this fix does not
   introduce a new fitted lambdaT3 or relax the scenario conditions.

## Install (Windows PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_General_Scalar_48_Campaign_Recovery_Fix.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_general_scalar_48_campaign.py
```

## Resume the actual calculation

```powershell
python -u -m scripts.RunT3GeneralScalar48Campaign `
  --run-id bulk48_v1 --resume --phase all `
  --sobol-points 8 --max-nfev 25
```

The completed per-case fit files and valid native beta exports are reused.
Missing or previously rejected cases (especially the E class) are attempted;
there is no intentional expensive re-execution of the successful A/B/C/D fits.
Do **not** delete the original `bulk48_v1` directory. The aggregate
`campaign.csv`, `campaign.json` and `summary.md` are refreshed.

Results:

```powershell
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\bulk_48\bulk48_v1\summary.md
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\bulk_48\bulk48_v1\errors.log -Tail 35
```

The aggregate now records `lambdaT3_input`, `c5_source_kind`,
`boundary_relaxed` and distinguishes total conditional accepted results from
those satisfying their original topology-coupling boundary.

## Scientific interpretation

This is a **conditional general-quartic numerical** comparison. It evolves
registered additional renormalisable quartics but omits permitted extra
holomorphic quadratic operators and complete general-scalar EFT1 Wilson
matching/running. The Matchete/RGBeta normalization is not physically
certified. Four extreme-hypercharge models use legacy one-scale C5 rather than
the hierarchical fixed-order C5 source. An accepted five-observable chi2 is
not a proof of predictivity, stability or a complete physical EFT.

In particular, B/C(-1) full-UV accepted results were obtained after changing
the historical lambdaT3 UV boundary; the summary labels them as modified
boundary fits. None of the other model/scenario constraints is relaxed.

## Recovery after 19/48 completed: C(+1) search bound and stale errors

The 19-accepted bulk summary is an actual outcome: 17 original-boundary
conditional candidates and two explicitly modified-boundary B(-1)/C(-1) UV
candidates. Do not interpret `NoMassSafeNeutrinoFit` as a failed oscillation
fit: such a candidate never reached the C5 evaluator.

C(+1) frozen-small was marked `FailedOrUnavailable` because its historical
Yukawa seed has an element with absolute value `2.746132227654911`, exceeding
the arbitrary `abs_yukawa_max=2.5` in the campaign's search box. The current
search cap is `sqrt(4*pi) = 3.5449...`, consistent with the already existing
maximum-Yukawa-singular-value diagnostic. This fixes a search-program error;
it does not waive the per-trial singular-value screen or guarantee acceptance.

`errors.log` remains an append-only audit trail. `errors_latest.log` is
replaced at the beginning of each new invocation and contains **only failures
from that invocation**. Old E-class errors therefore remain available for
audit but should not be treated as current failures.

### One targeted recovery command

This single-model command reuses saved native beta equations and 47 previously
reported case results. It does **not** rerun the other 15 models, nor change
the physical benchmark boundary.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -m pytest -q tests\python\test_t3_general_scalar_48_campaign.py
python -u -m scripts.RunT3GeneralScalar48Campaign `
  --run-id bulk48_v1 --resume --phase fit --models C_+1 `
  --sobol-points 16 --max-nfev 60
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\bulk_48\bulk48_v1\summary.md
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\bulk_48\bulk48_v1\errors_latest.log
```

The campaign always regenerates its complete 48-row `campaign.csv` summary,
even when a subset is selected. Existing `fit.json` records are reused and
not overwritten. Running the same command a second time will reuse the newly
written C(+1) result; do not expect an automatic second optimisation.

Further poor-fit or mass-screen cases require a separately specified numerical
search strategy or a clearly labelled modified-boundary study. They should
not be silently reclassified as accepted.
