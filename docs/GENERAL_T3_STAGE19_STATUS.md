# Stage 19 — selected extra-quartic Weinberg projection and low-energy transport

## Established inputs

- Stage 16: B/C(-1) grouped, fermion-first UV and scalar-only EFT RGBeta running, including all six special quartics, and **baseline** physical hierarchical C5.
- Stage 17: **new** quartic-dependent scalar-pair part of the one-loop dimension-five EFT1 Wilson-tensor anomalous dimension, based on Matchete PL LLSS tree seed.
- Stage 18: numerical integral of this **selected source**, using running EFT1 renormalisable couplings and *frozen* tree LLSS Wilson seed at mu_F. Four nonzero mixed-lepton Higgs-real components in both models (actual Windows runs).

## New verification and calculation

**Symbolic closure:** Use the project's authoritative `RGE.matching.WeinbergTensorAdapter.build_weinberg_wilson_tensor(kappa=1)` to generate the same SU(2) real-scalar template that validates the production EFT1 Wilson RGE. Independently recalculate the extra scalar-pair beta from the actual Matchete PL Wilson seed and the two SU(2) potential tensors for **all 30 independent lepton-doublet/Higgs-real index combinations** (`i<=j`, `a<=b`). Require exact identity

    delta_beta^(1)_{ijab} = beta_kappa^(1) * WeinbergTemplate_{ijab}

and require the four mixed-lepton components to agree with the saved Stage 17 expressions. If an index/sign/phase fails, **no numerical C5 is exported**. This check must be performed with the *actual local* seed; it has not been run here.

**Physical normalisation:** Independently compare the existing successful `eft1_wilson_rge.json` `weinberg_subspace_validation.beta_kappa_16pi2` with the `final_weinberg_coefficient.json` direct-running expression after inserting the canonical one-generation `C12 = conjugate(y1)*conjugate(y2)/MF` and extracting its coefficient of `log(MS/MF)`. This establishes that the RGE tensor's `kappa` uses the same one-loop direct-running C5 convention as the established physical matching adapter.

**Numerical projection:** Given the certified complex Weinberg template and the four integrated Stage 18 `C_{nu,e,a,b}`, the mixed-lepton template ratios are

    (a,b)=(1,3): -1/4
          (1,4): -i/4
          (2,3): -i/4
          (2,4): +1/4

for polynomial `kappa*(nu*H0 - e*H+)^2` and tensor `C=(1/4) fourth derivative`. Project by least squares and verify the *complete matrix* residual is smaller than 2e-8 relative, rather than using merely the four component norms.

**Scalar boundary:** Since the selected source is already LLHH with no S1 or S2 legs, the tree-level projection at the scalar threshold simply carries its coefficient to SM+Weinberg. **Finite one-loop scalar-threshold effects are not computed.**

**SM running:** Pass this **selected-source-only** symmetric 3x3 matrix through the existing `T3Trajectory.continue_with_weinberg_running`, from the grouped scalar matching scale `mu_S=1100 GeV` to `mu_low=246.22 GeV`, and convert each endpoint by `delta m_nu = -(246.22 GeV)^2/2 * delta C5`. Saved outputs remain separate from Stage 16 baseline.

## Deliberate limitations

- Full coupled EFT1 LLSS Wilson evolution and closure are not integrated (source-only, frozen LLSS seed).
- F and S matching remain **grouped**; individual heavy and scalar mass thresholds are not crossed.
- No finite extra contribution at the F or scalar threshold is included.
- With extra quartics generated from UV zero boundary, this contribution starts at **formal O(hbar^2)**; this is **not** a complete two-loop calculation, and the low-energy source-only output is **not a full physical neutrino mass prediction**.
- The running baseline and production matching results are not overwritten, nor are the source or GitHub files modified.

## Install and run (PowerShell, repository root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch19.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_stage19.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage19.ps1
```

Stage 19 uses existing Stage 17/18 user-generated reports and the *production* Matchete/RGBeta exports under `output/T3_B_alpha_m1/data/` and `output/T3_C_alpha_m1/data/`; no fresh symbolic matching or Stage 18 quadrature is required. It reruns two short renormalisable RGBeta numerical segments to recover the current consistent SM matching boundary.

## New files

- `Numerical/experimental/GeneralT3WeinbergProjection.py`
- `Numerical/experimental/GeneralT3ExtraLowEnergyTransport.py`
- `scripts/RunGeneralT3Stage19.py`
- `scripts/RunGeneralT3Stage19.ps1`
- `tests/python/test_general_t3_stage19.py`
- `docs/GENERAL_T3_STAGE19_STATUS.md`

## Output

- `output/general_t3_stage19/B_m1_selected_extra_weinberg_transport.json`
- `output/general_t3_stage19/C_m1_selected_extra_weinberg_transport.json`
- `output/general_t3_stage19/stage19_summary.json`

A failure writes an exact error to `stage19_summary.json` and exits nonzero, without manufacturing a `delta C5` value.
