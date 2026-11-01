# T3 general 16-model: large verification and reporting batch

## Scope

This adds a **read-only** post-processing route for the already completed general-interaction 16-model Matchete/RGBeta run. It does not edit `pipeline.py`, the generated native results, the model definitions, or GitHub.

It joins three tasks in one pass:

1. Re-evaluate the source-fingerprinted **structural inventory audit** (`AuditFull16GeneralSymbolicPhysics.audit_model`). The existing detailed Matchete checks and native RGBeta keys must still be consistent.
2. For **UV and EFT1**, independently derive and extract exact one-loop gauge factors, using native RGBeta `report_betas`. Tested structures include: gauge self-beta coefficients `b_Y,b_2,b_3`, diagonal scalar quartic and mass gauge terms `g_i^2 x`, and UV `y1,y2,y1Conj/y2Conj` matrix gauge terms. Additional quartics are handled by the declared field-count vectors, including the 12 previously registered RGBeta couplings used by B/C(-1). Exact arithmetic is via `Fraction`, SymPy and the Mathematica InputForm parser.
3. Generate the **physical coupling dictionary worklist** across 50 extra quartics, eight extra scalar mass bilinears and 38 Yukawa operators, without inventing signs/normalisations. Multi-invariant scalar monomials are labelled as undetermined basis blocks.

**This does not independently certify:** (a) the full Matchete↔RGBeta action-level map (relative Lagrangian/potential sign, factorials, complex conjugation, phases, native SU(2) component intertwiners), (b) all non-gauge beta-function coefficients, (c) extended EFT1 Wilson-operator running or general threshold matching / final C5. B(-1) and C(-1) still require mixed-scalar mass diagonalisation for matching.

## Files to place in the project

- `scripts/T3GeneralGroupFactorProjection.py` — independent exact group-factor formulas and Mathematica InputForm gauge-term extraction.
- `scripts/RunT3GeneralBatchVerification.py` — read-only orchestration, provenance checking, operator dictionary and report writing.
- `tests/python/test_t3_general_group_factor_projection.py` — six offline tests and negative controls.
- `docs/T3_GENERAL_GROUP_FACTOR_BIG_BATCH.md` — this instruction document.

These files are **new** and do not overwrite previous symbolic modules.

## Commands (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$env:PYTHONUTF8 = "1"
python -X utf8 -m pytest -q tests/python/test_t3_general_group_factor_projection.py
python -X utf8 -m scripts.RunT3GeneralBatchVerification --self-test

# Optional pilot, uses the native A(-2) files from the normal 16-model run:
python -X utf8 -m scripts.RunT3GeneralBatchVerification --model A_-2

# Full batch; no Wolfram calculation or numerical fitting:
python -X utf8 -m scripts.RunT3GeneralBatchVerification
```

Default input:

`output/studies/general_scalar_extended/complete_scalar_symbolic/full16_general_v1/`

Default output:

`Reports/output/full16_general_symbolic/GeneralBatchVerification/`

Files created:

- `README.md` — 16-model comparison with links into the existing full symbolic Lagrangian/RGE reports.
- `models/A_m2.md`, etc. — one model's formulas, exact gauge-factor comparison and operator worklist.
- `group_factor_coefficients.csv` — each exact coefficient comparison: model, stage, beta, gY/g2/g3, independent prediction, native projection, difference and status.
- `scalar_gauge_structure.csv` — independent SU(2) and U(1) scalar gauge wavefunction factors by interaction.
- `physical_coupling_dictionary.csv` — each extra scalar/mass/Yukawa interaction and its physical-map certification status.
- `native_beta_term_index.csv` — complete native UV/EFT1 beta-term index across the study (verbatim terms; non-gauge coefficients **not** independently certified).
- `model_status.csv`, `summary.json` — machine-readable checks and limitations.

The existing comprehensive general BSM Lagrangians and native UV/EFT1 beta expressions remain at:

`Reports/output/full16_general_symbolic/GeneralSymbolic/models/`

## Independent exact group-factor formulas

In conventions `Q=T3+Y`, `T(2)=1/2`, with three heavy fermion flavours:

\[
C_2(d)=\frac{d^2-1}{4}, \qquad T(d)=\frac{d(d^2-1)}{12}.
\]

For scalar monomial `O=\prod_i \phi_i`, its diagonal one-loop gauge-linear coefficient is

\[
\left.16\pi^2\beta_{\lambda_O}\right|_{\text{linear gauge}}=
-3\lambda_O\left[g_2^2\sum_i C_2(R_i)+g_Y^2\sum_i Y_i^2\right].
\]

This statement concerns the **diagonal** gauge-linear projection only; it does not determine non-gauge invariant-mixing and quartic recoupling coefficients.

For the three Yukawas' gauge term, the two fermion legs are `L` and `F`:

\[
\left.16\pi^2\beta_y\right|_\text{gauge}=
-3y\left[g_2^2\{C_2(L)+C_2(F)\}+g_Y^2\{Y_L^2+Y_F^2\}\right].
\]

The gauge-beta coefficients are calculated from the SM (`41/6`, `-19/6`, `-7`) plus the two complex BSM scalar multiplets and three massive F flavours. The F multiplet contributes one Weyl flavour per family for neutral self-conjugate representations, otherwise two for vector-like Dirac fermions.

As an exact **A(-2)** reference, `S1=(1,-1)`, `S2=(3,0)`, `F=(2,-1/2)`:

| Quantity | UV | EFT1 after F |
|---|---:|---:|
| coefficient of `gY^3` in `beta_gY` | `55/6` | `43/6` |
| coefficient of `g2^3` in `beta_g2` | `-1/2` | `-5/2` |
| coefficient of `g3^3` in `beta_g3` | `-7` | `-7` |

## How to interpret statuses

- `VerifiedRGBetaGaugeCoefficient`: exact numerical-rational agreement between an independent representation prediction and the projected saved RGBeta expression **for that specific gauge monomial**.
- `GaugeCoefficientMismatch`: exact values disagree. Inspect the corresponding `report_betas` original entry; this could reveal a source/model-convention discrepancy.
- `NotProjected`: the native term could not be parsed *unambiguously*. The test does not silently claim a pass. Investigate before asserting a complete group-factor check.
- `NotCertified`: physical tensor-basis action-level map and/or non-gauge coefficient equality is still open, **even if every gauge-factor row passes**.

A successful batch does not require regenerating the 16-model Matchete or RGBeta files. If any source fingerprint is stale, the structural gate fails closed.
