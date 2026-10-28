# T3 scalar relevance — Batch 4: one-loop quadratic RG source audit

## Goal and fixed scope

This batch follows the local 16/16 Wick-contraction audit and the 16-model
quadratic mass-spectrum study. Its specific question is whether the *allowed
additional quadratic mass coefficients* can consistently be set to zero when
the other scalar couplings are nonzero.

The working source is the existing `a02e04c6def5cf8522b2f612ccd6bff7d9963b91`
general-physics branch. No production Matchete or RGBeta source is modified.

Run an offline, **read-only** analysis of the 16 saved quadratic `UV` and 16
saved `EFT1` beta exports, plus 16 saved experimental EFT1 boundary/running
reports in `output/general_t3_quadratic_eft1/`. Unlike a live RGBeta run, the
audit merely analyzes the saved results. All 48 source files are SHA256
fingerprinted. Missing or incompatible inputs cause failure.

## Physics

Adopt the saved output convention

\[
16\pi^2\frac{d\mu_2}{d\ln\mu}=\beta^{(1)}_{\mu_2}.
\]

For example, the candidate native RGBeta mass-beta expression for B(-1) and
C(-1) has a source at zero mixed bilinear

\[
\left.\beta^{(1)}_{\mu_{12}^2}\right|_{\mu_{12}^2=0}
=6\lambda_{S_1S_1^{\dagger 2}S_2^\dagger}\,m_{S_1}^2
 +6\lambda_{S_1^\dagger S_2S_2^{\dagger 2}}\,m_{S_2}^2.
\]

These are respectively the saved RGBeta symbols
`lambdaS1S1bar2S2bar` and `lambdaS1barS2S2bar2`. Both are additional
symmetry-allowed quartics. Therefore **zero quadratic mixing alone is not a
generically RG-invariant restriction if these quartics are nonzero**.

Similarly, for A(+0), the candidate mass beta gives

\[
\left.\beta^{(1)}_{\mu_{S_1}^2}\right|_{\mu_{S_1}^2=0}
=6\lambda_{\mathrm{Q001300I1}}m_{S_1}^2
 +6\lambda_{\mathrm{Q000211I1}}m_{S_2}^2,
\]

with symbols `lambdaT3GeneralQ001300I1` and
`lambdaT3GeneralQ000211I1`. Other allowed bilinears have analogous terms.
The eight quadratic operators include **six neutral self-bilinears** and
**two B/C(-1) mixed-doublet bilinears**.

Importantly, when the additional quartic coefficients in each one-loop mass
beta are *also* set to zero, the explicit quadratic beta-function source
vanishes. This is a **mass-beta source statement only**; it does not prove
that the entire reduced quartic+quadratic system is RG closed.

Indeed, `output/general_rgbeta_source_proof.json` already shows that at least
two B/C(-1) Higgs-bearing extra quartics can be regenerated from the ordinary
`lambdaT3`, `y1`, and `y2` couplings. Starting from *all* such extra couplings
zero, a generated extra quartic may subsequently source a mass bilinear at an
additional loop order. That mechanism requires a full coupled assessment
before a globally RG-invariant ansatz can be claimed.

## Fixed-slope numerical diagnostic

The new script reads each candidate EFT1 initial state at `mu_F` and evaluates

\[
\Delta\mu_2^{\mathrm{LL,\ fixed\ slope}}
\simeq \frac{\ln(\mu_S/\mu_F)}{16\pi^2}
 \beta^{(1)}_{\mu_2}\Big|_{\mu_2=0,\;\mathbf g=\mathbf g(\mu_F)},
\]

with `mu_S < mu_F`. It records
`|Delta mu2|/min(mS1Sq,mS2Sq)` as a **dimensionless benchmark diagnostic**.
The source quartics and diagonal masses are held fixed at their saved EFT1
initial values. This is **not** a newly integrated zero-bilinear trajectory,
not a difference between two solutions with identical full UV boundary
conditions, and not a physical `C5` or neutrino-mass correction.

The saved quadratic EFT1 trajectories start from
`mu2_fraction=0.01`; their evolution **cannot** be used as a zero-boundary
solution without re-solving the coupled system. The script deliberately
calculates a local leading-log slope at `mu2=0` instead.

## How the source is verified

The script substitutes `mu2` and `Bar[mu2]` with zero in the complete saved
beta expression using a restricted arithmetic syntax tree. It checks that:

1. The remaining expression contains only the two diagonal scalar masses and
   allowed additional quartic symbols in that model.
2. Setting those additional quartics to zero annihilates the mass beta source
   algebraically, not merely to numerical tolerance.
3. UV and EFT1 source expressions agree after substitution in the saved
   output convention.
4. The EFT1 boundary has all required source coupling values, positive mass
   parameters and a descending positive scale interval.

Unsupported functions, novel source symbols, missing exports and inconsistent
metadata fail closed. No `eval` is used.

## Run (PowerShell)

Extract the ZIP into the **local** project root, preserving paths:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Scalar_Relevance_Batch4.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_quadratic_rg_closure_batch4.py
python -m scripts.AuditT3QuadraticRGClosure
Get-Content -Encoding UTF8 output\scalar_relevance\quadratic_rg_zero_boundary_summary.md
```

### New project files

- `scripts/AuditT3QuadraticRGClosure.py`
- `tests/python/test_t3_quadratic_rg_closure_batch4.py`
- `docs/T3_SCALAR_RG_CLOSURE_BATCH4.md`

### Generated locally by the script

- `output/scalar_relevance/quadratic_rg_zero_boundary_audit.json`
- `output/scalar_relevance/quadratic_rg_zero_boundary_audit.csv`
- `output/scalar_relevance/quadratic_rg_zero_boundary_summary.md`

## What remains before a final physical omission claim

The two main outstanding tasks are **coupled zero-boundary RG closure** for
all omitted quartics and bilinears, and checking the effect of nonzero allowed
scalar mass mixing on the actual one-loop `C5` in a physically normalized
mass basis. Omission of tree-level allowed bilinears cannot be justified
solely by loop-counting unless their zero-boundary value and radiative
regeneration are both controlled.

**Physical restrictions:** Matchete–RGBeta coupling and real-field mapping
remain uncertified; the complete 16-model mass-eigenstate matching is not
implemented. Nothing in this batch certifies a physical neutrino mass.
