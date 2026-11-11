# Worker C handover — 2026-11-12

## Files (new, additive)
- `Numerical/joint_scalar_matching/worker_c_joint.py`
- `tests/test_t3_worker_c_joint.py`

The package relies on existing `Numerical/joint_scalar_matching/core.py` and does not replace it. The code was authored against public GitHub `main` inspected on 2026-11-12, **not** the inaccessible Windows working copy. Check local `core.py` exposes `match_c5`, `check_tensors`, `EigenBasis` before installing.

## Verified public implementation inventory
- `core.py`: finite MSbar zero-momentum scalar-pair bubble `-C_LLSS*lambda_HHSS*B0bar/[2(16pi²)]`, with full real mass diagonalization, mixed-eigenstate loops, and strict Weinberg projection.
- `adapter.py`: canonical `(3,3,2,2,n,n)` complex `C_LLSS`, `(4,4,n,n)` real `lambda_HHSS`, `(n,n)` positive real scalar mass matrix, `(3,3)` inherited `C5_EFT1`; convention `T3_REAL_LLSS_HHSS_WEINBERG_v1`.
- `README.md`: other finite contributions beyond this bubble have **not** been certified; do not relabel a scalar bubble as full one-loop threshold matching.
- `README_PHYSICAL_READINESS_BATCH_V3.md`: signed Matchete GammaCC -> Weyl normalization was blocked in the documented native Feynman rule diagnostic; tensor provenance is mandatory.

Source URLs:
- https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/blob/main/Numerical/joint_scalar_matching/core.py
- https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/blob/main/Numerical/joint_scalar_matching/adapter.py
- https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/blob/main/Numerical/joint_scalar_matching/README.md
- https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/blob/main/Numerical/joint_scalar_matching/README_PHYSICAL_READINESS_BATCH_V3.md

## Interface

```python
from Numerical.joint_scalar_matching.worker_c_joint import match_joint_scalars
result = match_joint_scalars(eft1_wilsons, scalar_background, mu_s=1000.0)
```

`eft1_wilsons`: `C5_EFT1` already evolved to mu_S, `C_LLSS` already evolved to mu_S, plus `mu_s_gev`, `tensor_convention`, `renormalization_scheme`, and `provenance={'source': ...}`.

`scalar_background`: `mass_squared_gev2`, `lambda_HHSS`, `scalar_groups=['S1', ...,'S2',...]`, `scalar_mass_origin`, same metadata; optionally `lambda_HHSS_channels` with matching unique `quartic_channel_names`.

Returns `c5_inherited`, `c5_finite_scalar`, `c5_total`, pair and quartic-channel breakdown, the rotation, assumptions and omitted-term inventory. Returns NumPy arrays, not JSON; Worker D can serialize them using existing project encoders.

## Scope and deficiencies

**Computed:** all real components of both complex multiplets integrated jointly; one-insertion finite `LLSS x HHSS` bubble. Inherited C5 is added once without any manual EFT1 log.

**Not computed:** independent one-loop scalar topologies/other EFT1 operators, scalar external-leg matching or finite counterterms where applicable. No Matchete scalar export was transformed into a new finite coefficient; such a result would require validated amplitude/operator normalization and UV-minus-EFT subtraction. The existing Matchete script inventory includes `Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl`, `Lagrangian/RunGeneralDMinus2ScalarThreshold.wl` and `Lagrangian/RunGeneralDMinus2SequentialScalarThreshold.wl`; their complete generated output has not been verified against the current local repository. Therefore no claim of complete one-loop scalar matching.

## Installation (PowerShell)

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive -LiteralPath "$HOME\Downloads\T3_worker_C_joint_20261112.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p 'test_t3_worker_c_joint.py' -v
```

No historical outputs are modified and the test uses in-memory synthetic tensors. Before a physical run Worker D must insert Worker B's **same-scale normalized** mass and quartic tensors and the coupled solver's **same-scale evolved** Wilson arrays. Existing upstream EWSB-resummed mass convention, if used, must be labelled explicitly.
