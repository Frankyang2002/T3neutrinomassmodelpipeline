# Research General Scalar And Pilots

> Historical research record. This document combines **original source
> Markdown** under original-path headings. Local Markdown links to moved
> files are updated; exact original bytes are in the external archive.
> This is not a new certification
> of physics, a rewritten conclusion, or a new computational workflow.
> Original bytes and SHA-256 values are stored in the external archive.

## Contents

- [16 × 3 general-scalar campaign — result-recovery fix](#source-t3-general-scalar-48-campaign-md) — `docs/T3_GENERAL_SCALAR_48_CAMPAIGN.md`
- [General-scalar B/C pilot against the three historical studies](#source-t3-general-scalar-bc-historical-pilot-md) — `docs/T3_GENERAL_SCALAR_BC_HISTORICAL_PILOT.md`
- [Batch 13 — three-study general-scalar numerical campaign (B/C-supported scope)](#source-t3-general-scalar-numerical-campaign-md) — `docs/T3_GENERAL_SCALAR_NUMERICAL_CAMPAIGN.md`

---

<a id="source-t3-general-scalar-48-campaign-md"></a>

## Original document: `docs/T3_GENERAL_SCALAR_48_CAMPAIGN.md`

**Original SHA-256:** `337602699a91b19d1d8dc6548230bbb6b8d86a73ef2710e93a7bf6e9b669ea8e`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: 16 × 3 general-scalar campaign — result-recovery fix</summary>

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

</details>

---

<a id="source-t3-general-scalar-bc-historical-pilot-md"></a>

## Original document: `docs/T3_GENERAL_SCALAR_BC_HISTORICAL_PILOT.md`

**Original SHA-256:** `54cf455d79e716f5d8debfa005d355ab89be5ef31616f8b18d9a1312e99f7b76`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General-scalar B/C pilot against the three historical studies</summary>

# General-scalar B/C pilot against the three historical studies

## Scope

A **separate** exploratory comparison using original saved fitted B(-1)/C(-1)
Yukawa, fermion-mass and scalar inputs. The program reads each source under:

- `output/studies/lambda_5e-9_refit/{fitted,fits}/` — small \(\lambda_{T3}=5\times10^{-9}\), **UV RG on**;
- `output/studies/fermion_start_lambda_comparison/5e-9/{fitted,fits}/` — same topology quartic, **UV frozen**, boundary at fermion threshold;
- `output/studies/fermion_start_lambda_comparison/0p5/{fitted,fits}/` — \(\lambda_{T3}=0.5\), **UV frozen**, boundary at fermion threshold.

No inputs, fit files, historical output, Stage 36, or Batch 7–8 report is modified.
The generated study output is newly created only under
`output/studies/general_scalar_extended/<case>/<run_id>/`. An existing run ID
is rejected; use a distinct `--run-id` for any new scenario or repeated run.

## Benchmarks and actual calculations

For each of B(-1) and C(-1) the operator registry has **seven** additional
complex quartic slots (six existing RGBeta coordinates and one generated
sidecar) and one allowed mixed quadratic operator.

Use the original fitted ordinary scalar boundary, gauge/Yukawa/fermion values
unchanged. At the declared input boundary:

| Scenario | 7 additional quartics | mixed quadratic coefficient |
|---|---|---|
| `restricted_zero` | all exactly 0 | 0 |
| `general_quartics` | each `+0.1 abs(lambdaT3)` by default (real) | 0 |
| `general_quartics_and_mixing` (`--include-mixing`) | same nonzero quartics | `+0.01 min(mS1Sq,mS2Sq)` by default (real, GeV²) |

These are **chosen illustrative candidate RGBeta-coordinate initial values**,
not fitted values or established physical normalizations. These choices do not
impose positivity/vacuum stability on the full quartic potential.

For UV-on: integrate the **full coupled** candidate RGBeta system from 10^7 to
10^5 GeV, tree-project through fermion decoupling, and integrate the full
candidate EFT1 system from 10^5 to 10^3 GeV. For UV-frozen: do **not** run the
UV segment; impose the original fitted input values plus the chosen additional
couplings at 10^5 GeV, project fermions, and run only EFT1 to 10^3 GeV. No
fermion-threshold finite matching of scalar couplings is included.

Compare diagonal running scalar mass parameters, the topology quartic,
all seven additional quartics and (optionally) the scalar-mass eigenvalues and
mixing parameter from the specified SU(2) candidate mass matrix.

## One-scale C5 preview (optional)

`--c5-preview` reads the existing B/C one-scale Matchete ordered LLHH
coefficient and evaluates the symmetrized 3x3 C5 twice. **All** of the
Yukawas, fermion masses, scalar masses, scalar mixing matrix (mu2=0),
\(\lambda_{T3}\) and Matchete symbolic source are held fixed at the
*restricted* scenario values **in an explicitly hypothetical unmixed projection**. Only the two direct additional quartics
are varied to the `general_quartics` EFT1 endpoint values. The diagnostic is

\[
  \epsilon_{C_5}^{\rm identity\ preview} =
   \frac{\|C_5(\lambda_{a,b}^{\rm RGBeta})-C_5(0,0)\|_F}
        {\|C_5(0,0)\|_F}.
\]

**Crucial distinction:** even if the added quadratic coupling is **zero at the
input boundary** of `general_quartics`, coupled one-loop RGEs can induce a
**nonzero quadratic mixing at the EFT1 endpoint**. This is not a failed RG
run. The `--c5-preview` calculation explicitly sets mu2=0 in *both* C5
evaluations to isolate direct quartic insertion at a common unmixed background.
The actual endpoint quadratic coefficients and |mu2|/min(m²) ratios for both
scenarios are reported in the per-model JSON and the summary. No physical C5
for the mixed-mass endpoint is computed; therefore the preview must not be
interpreted as the total C5 difference including mixing.

**Other limitations:** this assumes the **unverified** identity mapping of the
nonzero RGBeta quartics to the Matchete symbols, freezes original study
Yukawas (not UV-evolved Yukawas), and substitutes a simultaneous one-scale
coefficient rather than performing sequential EFT1 Wilson matching. This is
**not** a physical \(C_5\) error estimate, neutrino mass prediction, or
complete general-scalar matching. The *actual* nonzero quadratic mixing is not inserted into the one-scale
\(C_5\) evaluation: it can affect leading one-loop matching and needs a
certified mass-basis/field-map treatment.

The historical fit chi² and `Success`/`PoorFit` flags are only provenance;
new scalar scenarios are not refitted. Do not label a poor historical fit as
physically accepted because the numerical RG solver succeeds.

## Installation and PowerShell execution

Extract the ZIP at the **project root**, keeping the three historical study
folders unchanged. The archive contains **full new files**, no code patches.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

# Test without rerunning the fitted studies:
python -m pytest -q tests\python\test_t3_general_scalar_bc_pilot.py

# Fast pilot: B(-1) with the small UV-frozen study only:
python -m scripts.RunT3GeneralScalarBCPilot --case small_uv_frozen_5e-9 --model B_-1 --run-id pilot_q0p1

# All three historical cases, both models; baseline versus 7 extra quartics:
python -m scripts.RunT3GeneralScalarBCPilot --all --c5-preview --run-id pilot_q0p1_c5_unmixed_projection

# OPTIONAL: mixed-mass scenario also, as a separate NEW result ID:
python -m scripts.RunT3GeneralScalarBCPilot --all --include-mixing --run-id pilot_q0p1_mu0p01

# Read the six summaries:
Get-ChildItem .\output\studies\general_scalar_extended\*\pilot_q0p1_c5\summary.md | ForEach-Object {
    Get-Content -Encoding UTF8 $_.FullName
}
```

Default extra-quartic fraction is 0.1 of the original fitted |lambdaT3|;
`--quartic-fraction 0.01` selects 1%. `--mu2-fraction 0.01` sets 1% of the
smaller initial diagonal mass squared. Both use **positive real** coefficients.
For distinct choices use distinct `--run-id` values.

## Verification status

Offline tests exercise input/scenario validation, path exclusivity, candidate
mass spectra, RG endpoint comparisons and C5 preview guards. Full numerical
integration on the user's Windows checkout is required before quoting any
actual pilot values. The implementation reuses the saved candidate RGBeta
beta exports; no Wolfram regeneration or GitHub changes are necessary.

## Windows zero-input-mixing preview correction

The first runner incorrectly required `quadratic_mu2_gev2 == 0` at the **EFT1
endpoint** to run the optional `--c5-preview`. A vanishing boundary coefficient
can become nonzero under the coupled RGEs. The corrected runner does not erase
or falsify the generated coefficient: it records it in each RG scenario and in
`C5_preview.unmixed_projection_info`, while the *separate direct-quartic*
preview explicitly evaluates a hypothetical diagonal scalar mass matrix.
This does not provide a new full-mixing C5 calculation.

</details>

---

<a id="source-t3-general-scalar-numerical-campaign-md"></a>

## Original document: `docs/T3_GENERAL_SCALAR_NUMERICAL_CAMPAIGN.md`

**Original SHA-256:** `8ba038f575e8408e28493c43f42d04cde500eb5111749c43b44e52b59dab0dbd`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Batch 13 — three-study general-scalar numerical campaign (B/C-supported scope)</summary>

# Batch 13 — three-study general-scalar numerical campaign (B/C-supported scope)

## What the current code genuinely supports

There are **three** established historical studies and **16** ordinary T3 assignments: **48 model–study cells**. The complete supported candidate general-scalar oscillation-fit path currently exists for only **B(-1), C(-1)**: **6 cells**. The other **42 cells** (14 models per study) are deliberately marked `BlockedGeneralMatching`.

**Do not report the 42 blocked entries as new general-scalar fits.** Model-dependent full general-potential C5 expressions and appropriate quadratic operators are not incorporated for them. Existing 16-model Stage 31/34/35/36 calculations use earlier truncated/intermediate matching assumptions and do not remove that limitation.

The supported B/C path uses:

- Native candidate RGBeta UV/EFT1 with all seven B/C additional quartics and the holomorphic quadratic mass coefficient.
- The exact historical UV-on versus UV-frozen boundary choice.
- Matchete's **saved one-scale** mixed-scalar ordered C5 expression and Stage 23/26 **conditional** direct-quartic sign `lambda_M=-lambda_RGBeta`.
- Heavy fermion masses/Yukawas in the validated mass basis at `muF=100 TeV`.
- SM+Weinberg running from `muS=1 TeV` to `100 GeV`; existing five-observable NuFIT chi-square; existing 18-real-Yukawa optimiser.
- Independent illustrative scalar boundary: `quartic_fraction=0.1` times the absolute topology quartic; quadratic input `mu2_fraction=0` by default (nonzero values can be requested separately).

This is a *conditional exploratory* fit, not physical general-potential matching or a fully renormalised sequential EFT prediction. Fitting five observables with 18 free Yukawa parameters is not an independent prediction.

## Files (all new, full files)

- `scripts/RunT3GeneralScalarNumericalCampaign.py` — read-only full inventory plus optional candidate B/C fitting campaign.
- `tests/python/test_t3_general_scalar_numerical_campaign.py` — 22 offline tests for inventory, continuation, provenance and nonoverwrite policy.
- `docs/T3_GENERAL_SCALAR_NUMERICAL_CAMPAIGN.md` — this document.

No existing source file is replaced, no GitHub write is performed, and no historical study or prior Batch 9–12 output directory is modified.

## Commands (PowerShell, project root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_General_Scalar_Numerical_Campaign_Batch13.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_general_scalar_numerical_campaign.py

# Inspect 48 cells and detect accepted prior general-scalar B/C fits; WRITE NOTHING.
python -m scripts.RunT3GeneralScalarNumericalCampaign --plan

# Reuse any accepted and source-hash-matched B/C fits; refit only the missing/poor ones.
python -m scripts.RunT3GeneralScalarNumericalCampaign --run --campaign general_v1

# Optional: genuinely rerun all SIX supported B/C cells as NEW fit directories;
# computationally more expensive, never overwrites the results from general_v1.
python -m scripts.RunT3GeneralScalarNumericalCampaign --run --campaign general_fresh_v1 --fresh
```

The campaign uses existing Batch 11 fitter and Batch 12 18-parameter continuation when needed. If a first fit has chi-square above 5.86, continuation starts from the first fit's saved best Yukawas. A continuation result still above 5.86 is correctly recorded as a **poor exploratory fit**; there is no guarantee of acceptance.

Separate the `general_v1` and `general_fresh_v1` campaigns. Running either a second time with the same name **fails closed**. Each generated fit also has a new run ID and fails on destination collision. The accepted prior fit is reused only when the saved study/model, scalar boundary, certification flags, fitted configuration metadata, and **source file SHA-256 hashes** agree.

The full campaign report is:

```text
output/studies/general_scalar_extended/campaigns/<campaign>/campaign.json
output/studies/general_scalar_extended/campaigns/<campaign>/summary.md
```

Read it:

```powershell
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\campaigns\general_v1\summary.md
```

New B/C individual fits (if required) are stored at:

```text
output/studies/general_scalar_extended/<study>/<campaign>_B_m1_fit/
output/studies/general_scalar_extended/<study>/<campaign>_B_m1_cont/
output/studies/general_scalar_extended/<study>/<campaign>_C_m1_fit/
output/studies/general_scalar_extended/<study>/<campaign>_C_m1_cont/
```

The `*_cont` run only exists if the first fitting stage did not reach the chi-square cut. Each completed B/C fit writes its full fitted configuration, fit JSON and summary.

## Explicit block to completing all 48 general-scalar refits

For the other 14 assignments per study, the 38 new scalar-quartic slots and additional quadratic structures require appropriate model-specific EFT transfer and C5 response. The four extreme-hypercharge cases already use legacy C5 artifacts instead of the same hierarchical scheme in historical Stage 36. Neither the generic quartic RG state nor old model-specific C5 alone is a fully general-potential matching implementation. Implement, check, and document those C5 and quadratic interfaces *before* launching the missing 42 refits.

We intentionally stop rather than manufacture numerically plausible but incorrectly labelled general-potential results.

</details>

---
