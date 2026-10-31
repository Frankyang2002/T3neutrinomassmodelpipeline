# Workstream D — Stage 4 numerical validation and native state adapters

**Status:** Additive final numerical-interface tests. **Not** a new physical all-16 matching/RGE solution. No existing B/C runner, benchmark, output, symbolic model, or RGBeta beta export has been modified.

## Files / project-relative paths

1. `Numerical/experimental/GeneralT3NumericalStateAdapters.py` — native state <-> Stage 3 state-envelope round trip (UV, scalar-only EFT1, general RGBeta extra quartics, SM+Weinberg).
2. `Numerical/experimental/GeneralT3NumericalChecks.py` — conditional Wilson tensor (flavour and scalar axes), full complex symmetric C5, native Takagi residual, and B/C saved-report regression.
3. `scripts/AuditGeneralT3NumericalCompletion.py` — no-write-to-input dependency and B/C report validation; emits a separate summary and refuses overwrites.
4. `tests/test_general_t3_numerical_stage4.py` — synthetic structural, invalid-input, complete native round-trip and fail-closed tests, including an optional actual general-RGBeta registry test.
5. `docs/WorkstreamD_Stage4NumericalValidation.md` — this contract.

## Source interfaces inspected (read-only public repository)

- `Numerical.core.State` (`T3UVState`, `SMNumericalState`)
- `Numerical.core.StateVector` (`pack_uv_state`, `unpack_uv_state`, real/imag layout)
- `Numerical.running.IntermediateScalarState` / `IntermediateScalarStateVector` / `FermionThresholdBoundary`
- `Numerical.experimental.GeneralT3ScalarRGState` (`pack_general_scalar_state`, `unpack_general_scalar_state`, extra RGBeta quartics)
- `Numerical.running.SMWeinbergEvolution` (`SMWeinbergInitialConditions`, `K`)
- `physics.NeutrinoMass` (`numerical_neutrino_mass_matrix`)
- `physics.NeutrinoObservables` (`takagi_factorization`)
- `Numerical.experimental.T3BCCoupledWilsonTensor` (six independent lepton-flavour pairs)
- `scripts.RunT3BCEFT2Sequential` (saved conditional B/C report schema)
- Stage 3 `GeneralT3EFTHandoffController` and candidate readiness report.

The adapter calls the original production pack and unpack functions instead of recomputing RGBeta beta functions or choosing an independent field basis. For each native block, coordinates preserve the source naming, shape, dtype, and real/imag storage. Scales remain separate GeV metadata. Direct UV and EFT1 native states, general-extra-quartic native states, and final `SMWeinbergInitialConditions` are each supported.

For a general bilinear model, `mu2` is a distinct required complex scalar in GeV², reproducing the continuous solver's additional two-real-coordinate sidecar. This does **not** apply a quadratic mass-basis rotation; it only preserves its value.

A `NativeBridge.restore(...)` validates the **same-stage** numerical contract and delegates reconstruction to the original source unpacker. It does not support or imply matching into a different EFT. A matching callback remains mandatory for particle removal.

## Wilson / neutrino mass conventions

The B/C conditional Wilson tensor is `(3,3,N_s,N_s,6)` with `N_s=12` (both original scalars) or `N_s=8` (H and X), and must be symmetric under the two lepton indices and separately the two scalar indices. Its six-pair ordering is the project-defined `T3BCCoupledWilsonTensor.PAIRS`, not an independent three-generation gauge-index RGE. All Wilson coefficients require *explicit* operator-basis and normalization labels; providing labels is not physical certification.

The final SM+Weinberg source stores its complex symmetric Weinberg coefficient as `K`. The native adapter maps it to the controller `C5` block **without changing components**. The existing project convention is

\[
m_\nu=-\frac{v^2}{2}C_5,\quad v=246.22\ \mathrm{GeV},\quad [C_5]=\mathrm{GeV}^{-1}.
\]

`check_neutrino_matrix` calls **the existing** `numerical_neutrino_mass_matrix` and **the existing** `takagi_factorization` functions. It checks the Takagi residual and SVD masses independently. A separate C5 symmetry check is relative to the tiny matrix magnitude; the generic Stage 3 `1e-11` absolute symmetry tolerance alone is too permissive for realistic `~1e-15 GeV^-1` C5.

The saved B/C `mnu_at_100gev_ev` is in the charged-lepton mass basis, whereas `C5_at_100gev_gev_inverse` may not be. Consequently, the matrix elements need not match directly. The singular masses must match because charged-lepton basis rotations are unitary congruences; this is what the Stage 4 report check verifies.

## Regression values (conditional; saved source benchmark)

- B(-1): `[0.03853362292521312, 0.01561382145082619, 0.00045848474457153946]` eV
- C(-1): `[0.046405741900074446, 0.045923847867146636, 0.021684644219742994]` eV

These values, the Y/X ordering, finite-matching **missing** flags, and physical-certification **false** flags are checked against the existing saved or freshly rerun JSON outputs. Stage 4 never reconstructs missing matching coefficients or reinterprets them as numerical zero.

## PowerShell instructions

From the project root:

```powershell
python -m unittest discover -s tests -p "test_general_t3_numerical_stage4.py" -v
python -m unittest discover -s tests -p "test_general_t3_eft_handoff_controller.py" -v
python -m unittest discover -s tests -p "test_general_t3_bc_conditional_regression.py" -v
python -m scripts.AuditGeneralT3NumericalCompletion
```

The audit reads:

- `output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/workstream_d_fermion_audit_v2/audit_summary.json`;
- `.../workstream_d_handoff_readiness_v1/summary.json`;
- `.../bc_eft2_d_regression_v1/B_m1_eft2_sequential.json` and `C_m1_eft2_sequential.json`.

The audit writes only:

`output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/workstream_d_numerical_validation_v1/summary.json`

The script **requires** both real B/C report files rather than silently skipping absent data and does not overwrite an existing output. Use `--bc-reports` to select the historical `bc_eft2_sequential_v1` directory, or `--output` to choose a new report location for reruns.

If you want a fresh *independent* numerical regression from the original B/C runner, use an unused `--output-id` with `scripts.RunT3BCEFT2Sequential` and pass that output directory through `--bc-reports`. Stage 4 itself never reruns or modifies the original B/C pipeline.

## Workstream handoff status / limitations

The uploaded Stage 2/3 reports contain 84 candidate roots across all 16 models; B(-1) `fermion:0` and D(-2) `scalar:S2_minus` remain unbracketed. They are *candidate full-UV/grouped-EFT1 roots*, not physically recomputed below individual mass thresholds. Stage 3 has **zero registered RG and matching callbacks**.

After local Stage 4 passes, **Workstream D numerical infrastructure can be closed**. This is deliberately narrower than complete physically certified T3 predictions. Remaining physics is external:

- Workstream B: model-dependent EFT RGEs and beta maps after each removal.
- Workstream C: verified fermion/scalar matching (including finite B/C Y/X).
- Workstreams A/B/C: certified operator and scalar-field normalizations where still unresolved.

Those inputs must be provided before calling the universal controller for a fully physical sequential T3 spectrum.
