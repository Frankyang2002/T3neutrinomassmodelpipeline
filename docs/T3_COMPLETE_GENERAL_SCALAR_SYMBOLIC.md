# General-scalar T3: opt-in symbolic bilinear integration

<!-- T3_CURRENT_SCOPE_BEGIN: 2026-11-10 -->

**Document scope:** General symbolic exports exist; general numerical two-threshold Wilson matching is not implied.

**Current implementation status (10 November 2026):** The requested A(-4)/D(-2) **full-general, complete-one-loop, fermion-first -> joint-scalar -> SMEFT Yukawa fit is not yet an executable end-to-end numerical workflow**. The repository already implements general Matchete/RGBeta components and an actual symbolic `RunSequentialT3Matching` stage engine, but its general Wilson matching outputs are not fully connected to the two-model numerical fit. The current thesis Sobol pilot is **restricted/conditional** (analytic fermion triangle and selected scalar bubble), even if it reaches a low chi-square. See [implementation status](IMPLEMENTATION_STATUS.md), [A(-4)/D(-2) feasibility](T3_A4_D2_GENERAL_PIPELINE_FEASIBILITY_2026_11_10.md) and [implementation batches](T3_A4_D2_FULL_GENERAL_IMPLEMENTATION_BATCHES.md). Earlier numerical results and verification claims below refer to their **dated, historical configurations** unless an individual run manifest establishes otherwise.

<!-- T3_CURRENT_SCOPE_END -->
## Scope and scientific interpretation

This package adds the symmetry-allowed holomorphic scalar *quadratic* operators
into the **opt-in symbolic UV Lagrangian builder**. When enabled, the existing
independent quartic basis is included as well. The separate Matchete detailed
operator audit and RGBeta one-loop **UV and EFT1** beta generator are connected
by a new driver. All original restricted benchmark input files, outputs and
numerical physics implementations remain unchanged.

This is **not yet** a complete physical general-scalar sequential EFT calculation:

- B(-1) and C(-1) contain mixed scalar mass bilinears. Matchete's
  `HeavyMassBasis` diagnostic does not pass until the scalar fields and all
  couplings have been transformed to a consistent mass basis. The complete
  interaction expression is exported, but `RunModel.wl` is not allowed to
  continue to matching for those models.
- The diagonal/self-bilinear cases may assemble with `CheckLagrangian` passing;
  this does not certify a physically normalised mass-eigenstate `C5` matching
  calculation. No new `C5` is claimed by this package.
- The candidate RGBeta UV/EFT1 beta expressions are derived by the existing
  `RunGeneralT3QuadraticRGBeta.wl`, including both diagonal scalar mass betas,
  the new holomorphic mass beta, and full general quartic projections. The
  **Matchete/RGBeta field and coupling normalisation** remains uncertified.
- Full Wilson-operator anomalous dimensions in the intermediate EFT and final
  mass-eigenstate threshold matching are **not** implemented here. Do not use
  these exports as replacements for current fitted-neutrino results.

## Model-dependent additional quadratic mass operators

| Model | Bilinear type | Count |
|---|---|---:|
| A(0), D(0), E(0) | neutral S1 self-bilinear | 1 each |
| A(-2), D(-2), E(-2) | neutral S2 self-bilinear | 1 each |
| B(-1), C(-1) | S1/S2 mixed bilinear | 1 each |
| Other eight ordinary T3 models | none | 0 |

Every bilinear is represented by an independent **complex mass-squared**
coefficient plus its Hermitian conjugate. In the candidate RGBeta potential
convention, `V2 = mu2 * invariant + h.c.` without an extra factor of `1/2`.
The other 50 **model-slots** of extra quartics come from the existing
`GeneralScalarOperators.wl`/RGBeta registry; this package does not redefine them.

## Files and destinations

- `Lagrangian/interactions/T3Topology.wl` — **full replacement**, opt-in extra
  quadratic and quartic candidates.
- `Lagrangian/model/LagrangianBuilder.wl` — **full replacement**, fail-closed
  handling: do not discard a gauge-allowed but mass-basis-blocked bilinear.
- `Lagrangian/RunFullGeneralScalarBuild.wl` — **new** standalone Lagrangian
  export (not a matching computation).
- `scripts/RunT3CompleteGeneralScalarSymbolic.py` — **new** symbolic driver.
- `tests/test_t3_complete_general_scalar_symbolic.py` — **new** pure-Python tests.
- `docs/T3_COMPLETE_GENERAL_SCALAR_SYMBOLIC.md` — these instructions.

The existing `RunGeneralQuadraticMatcheteAudit.wl`,
`RunGeneralT3QuadraticRGBeta.wl`, `GeneralScalarOperators.wl` and
`general_rgbeta_quartic_registry.json` are required without modification.

## Installation and commands (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
# Back up the two replaced files before extracting:
Copy-Item Lagrangian\interactions\T3Topology.wl Lagrangian\interactions\T3Topology.wl.before_bilinears
Copy-Item Lagrangian\model\LagrangianBuilder.wl Lagrangian\model\LagrangianBuilder.wl.before_bilinears
Expand-Archive -LiteralPath 'C:\path\to\T3_full_general_scalar_symbolic.zip' -DestinationPath . -Force

python -X utf8 -m pytest -q tests/test_t3_complete_general_scalar_symbolic.py
python -X utf8 -m scripts.RunT3CompleteGeneralScalarSymbolic
python -u -X utf8 -m scripts.RunT3CompleteGeneralScalarSymbolic --model A_-2 --run --run-id a_m2_v1
python -u -X utf8 -m scripts.RunT3CompleteGeneralScalarSymbolic --all --run --run-id all16_v1
```

**Requirements:** working Mathematica/Wolfram FrontEnd + Matchete and RGBeta,
as required by the existing symbolic repository. The first two commands do not
invoke Wolfram. The `--run` calls do.

Outputs are written to:
`output/studies/general_scalar_extended/complete_scalar_symbolic/<run-id>/`.
Each model has a standalone Lagrangian interaction InputForm, detailed Matchete
audit and two RGBeta beta-function JSONs plus logs. `summary.json` and
`summary.csv` record status and failures. Existing run IDs are not overwritten.

For the **existing RunModel.wl** (ordinary matching driver), the flag is
`T3_EXPERIMENTAL_GENERAL_QUADRATICS=1`; leave it unset for all historical
numerical campaigns. Enabling this option without changing subsequent matching
**does not** certify the output as the fully general scalar neutrino coefficient.

## What remains before numerical use

1. Demonstrate field-basis and coupling-normalisation equivalence between
   Matchete and the RGBeta quartic/mass conventions.
2. For B(-1)/C(-1), diagonalise the mixed complex scalar mass matrix and
   transform all couplings before heavy scalar matching.
3. Derive/validate the general-bilinear one-loop `C5` for the relevant models;
   retain the old diagonal result only as a control.
4. Connect live generated UV/EFT1 scalar beta functions and matching into
   the numerical state and threshold evolution, then rerun benchmarks.

The first numerical re-fit or running comparison should be done only *after*
steps 1–3, to avoid giving the new bilinear incorrect physical effects.
