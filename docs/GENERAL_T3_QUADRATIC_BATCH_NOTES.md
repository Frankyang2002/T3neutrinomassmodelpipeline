# T3 general-scalar batch: eight quadratic operators and one-loop mass beta functions

## Scope / source checkpoint

Read-only source reference: `Frankyang2002/T3neutrinomassmodelpipeline`, public `main`, checked 2026-10-27, including:

- `docs/general_scalar_trial_operator_manifest.json`
- `Lagrangian/model/T3Fields.wl`, `Lagrangian/model/T3ModelCatalog.wl`
- `Lagrangian/interactions/GeneralScalarOperators.wl`, `GeneralScalarVerifiedQuartics.wl`
- `RGE/running/rgbeta/T3RGBetaModel.wl`, `RunGeneralT3RGBeta.wl`
- `Numerical/experimental/GeneralT3ScalarRGState.py`, `GeneralT3ScalarQuadraticMasses.py`
- `scripts/FermionThresholdHandoff.py`, `IntermediateEFTRunning.py`
- 16-model Stage29 registry and the Stage22–29 physical tensor-map documentation.

The user separately reports a successful 16-model **full quartic EFT1** computation. This batch does not rerun/rebuild the original matching or numerical infrastructure, and it does not overwrite Stage34–36.

## Verified independent operator inventory

For `Y(S1)=alpha/2`, `Y(S2)=(alpha+2)/2`, with complex Z2-odd scalars, the only additional dimension-two gauge singlets in the 16 ordinary models are:

| Model | Holomorphic bilinear, plus H.c. | Candidate mass symbol | Basis interpretation |
|---|---|---|---|
| `A_+0` | `(S1†)^2` | `muT3GeneralS1Sq` | singlet split |
| `A_-2` | `S2† · S2†` | `muT3GeneralS2Sq` | triplet split |
| `B_-1` | `eps(S1†,S2†)` | `muT3GeneralS1S2` | mixed doublets |
| `C_-1` | `eps(S1†,S2†)` | `muT3GeneralS1S2` | mixed doublets |
| `D_+0` | `S1† · S1†` | `muT3GeneralS1Sq` | triplet split |
| `D_-2` | `(S2†)^2` | `muT3GeneralS2Sq` | singlet split |
| `E_+0` | `S1† · S1†` | `muT3GeneralS1Sq` | triplet split |
| `E_-2` | `S2† · S2†` | `muT3GeneralS2Sq` | triplet split |

No extra bilinear appears in the other eight assignments.  **There are 8 model-specific complex slots, not 8 universal couplings in each model.**  In the real gauge-unbroken scalar component space, the bilinear contributes a symmetric holomorphic mass matrix `B`. The two B/C terms may be diagonalized after rewriting the conjugate scalar to the same electroweak representation. That is only a *basis choice*: Yukawa/quartic tensors must be transformed too, so omitting the interactions from the original declared field basis is not justified. The other six coefficients can have a removable phase, but their magnitude produces a physical splitting. No claim of a full theory basis-redundancy proof is made.

## What this batch actually implements

1. **Independent quadratic gauge-invariant inventory**: `scripts/AuditGeneralT3QuadraticBasis.py`, checks the original manifest without using it to generate the candidate list. Fail-closed 8/16 test.
2. **Matchete audit, experimental**: `Lagrangian/RunGeneralQuadraticMatcheteAudit.wl`, builds the unchanged ordinary Lagrangian, all experimental quartics and the one allowed quadratic (if any), and validates individually and collectively using `CheckLagrangian`. `scripts/RunGeneralT3QuadraticMatcheteAudit.py` orchestrates 16 models. Requires live Matchete front end.
3. **RGBeta one-loop mass betas, experimental**: `RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl`, adds each model's quadratic directly via the documented `AddScalarMass[..., SelfConjugate -> False]` API, alongside all Stage29 general quartics; obtains `BetaTerm` for the **new complex quadratic coefficient and the two diagonal masses**, plus the full projected quartic beta class, for `UV` or `EFT1`. The script does not modify baseline models.
4. **EFT1 coupled numerical mass/quartic RG flow, experimental**: `scripts/RunGeneralT3QuadraticEFT1.py` reconstructs the hashed Stage33 threshold boundary and Stage29 sidecars; evolves `mS1Sq`, `mS2Sq`, 38 new / 12 reused quartic slots across models, and **the 8 new complex mass coefficients** wherever allowed, subject to successful live beta export. It rejects missing beta keys and unsupported InputForm expressions. Optional `--also-export-uv` exports UV mass betas but does not integrate a generalized UV mass boundary.
5. **Full replacement Hessian library**: `Numerical/experimental/GeneralT3ScalarQuadraticMasses.py`, preserves the existing function signatures and uses independent analytic spectra and the direct complex potential as checks. The chosen spherical SU(2) basis, including the neutral triplet bilinear metric, is explicit.
6. **Offline tests**: `tests/python/test_general_t3_quadratic_full_batch.py` tests 16 zero-limit models, all 8 complex quadratic sectors with three complex phases, spectra, Hessian contractions, tachyon flags and manifest checks. No Wolfram numerical result is claimed without a local live run.

## Physical convention and important limitation

The RGBeta mass term uses the **potential** convention (RGBeta's Lagrangian scalar potential enters with a minus sign):

`V_extra = mu2 * (Sa†)^T C (Sb†) + h.c.`

No hidden `1/2` is inserted. In this convention, self-conjugate-direction branches are `m² ± 2|mu2|`, and the B/C two-doublet squared-mass branches are

`(m1²+m2²)/2 ± sqrt(((m1²-m2²)/2)² + |mu2|²)`.

**This does not certify the Matchete Lagrangian source sign or its triplet/Bar component convention.** The selected spherical SU(2) tensors are mathematical conventions. The Stage22–29 quartic mapping remains conditional. Matchete mass-parameter basis, mass eigenstate threshold matching, and the Wilson coefficients for the new interactions remain unresolved. RGBeta's API confirms `AddScalarMass` can declare a complex mass coupling, but its live execution and one-loop output must still be verified on the user's Wolfram installation.

## Installation

Extract ZIP preserving paths into:

`C:\thesistemp\T3neutrinomassmodelpipeline`

No existing production file is replaced. The *existing experimental* `Numerical/experimental/GeneralT3ScalarQuadraticMasses.py` **is intentionally replaced in full** (same exported function names, additional checks).

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Quadratic_Generalisation_Batch.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_quadratic_full_batch.py
python -m scripts.AuditGeneralT3QuadraticBasis
```

First *live RGBeta mass-beta + numerical EFT1* test:

```powershell
python -m scripts.RunGeneralT3QuadraticEFT1 --model B_-1 --mu2-fraction 0.01
```

Then all 16, optional UV mass beta exports:

```powershell
python -m scripts.RunGeneralT3QuadraticEFT1 --all --mu2-fraction 0.01 --also-export-uv
```

Independent *live Matchete* candidate-Lagrangian test:

```powershell
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --model B_-1
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --all
```

Output:

- `output/general_t3_quadratic_eft1/quadratic_basis_audit.json`
- `output/general_t3_quadratic_eft1/B_m1_eft1_quadratic_live_rgbeta.json`
- `output/general_t3_quadratic_eft1/B_m1_eft1_quadratic_running.json`
- `output/general_t3_quadratic_eft1/B_m1_uv_quadratic_live_rgbeta.json` (only with `--also-export-uv`)
- `output/general_t3_quadratic_eft1/all_models_quadratic_eft1_summary.json` (only after successful `--all`)
- `output/general_t3_quadratic_eft1/matchete/...` (only with Matchete audit)

The default `--mu2-fraction 0` is a **zero-extra-quadratic control**. With `--mu2-fraction 0.01`, the starting extra coupling in a model that permits one is `0.01 * min(mS1Sq,mS2Sq) * exp(i*phase)`; the other eight models have no such parameter. `--mu2-phase` is in radians.

## Verification gates and next physics work

- OFFLINE: syntax, 16-model operator count and Hessian spectra/direct complex potential checks.
- LIVE: Matchete `CheckLagrangian` must accept all eight alongside the 50 candidate quartic slots; RGBeta must return complete UV/EFT1 beta inventories for all 16 models; the extended solver must succeed across all 16 without missing keys.
- UNRESOLVED: independent Matchete/RGBeta physical real-field intertwiner, exact mixed-mass/barred-triplet normalization, generalized UV mass RG boundary, diagonalization at the scalar threshold, and new one-loop matching/Wilson mixing.

**A numerically successful EFT1 mass/quartic run is NOT proof of a physically complete general T3 neutrino mass model.** Production Stage34–36 results are not overwritten or reinterpreted.
