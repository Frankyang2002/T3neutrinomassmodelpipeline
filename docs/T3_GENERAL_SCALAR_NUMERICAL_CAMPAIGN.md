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
