# Stage 20 RGBeta v1.2.0 compatibility and candidate-gauge audit

## Failure identified

The local Wolfram log reported:

    General::invalidopt: Option CheckInvariance for function AddQuartic received invalid value True.
    General scalar quartic registration failed

The 16 generated `AddQuartic` expressions previously set `CheckInvariance ->
(Environment["T3_GENERAL_RGBETA_VERIFY_GAUGE"] === "1")`. When the user selected
`-CheckGauge`, that option evaluated to `True` and RGBeta v1.2.0 rejected the
registration. The Wolfram process therefore failed **before** computing betas.

## Repair

- The generator no longer passes the unsupported `CheckInvariance` option.
  The option uses RGBeta's default; it is not silently re-enabled.
- The `-CheckGauge` flag now performs an independent, explicit **structural**
  audit of all 38 *new* candidate quartic operators: U(1)Y neutral charge,
  Z2 parity, conjugated field orientation and SU(2) candidate multiplicity.
  The audit **is not RGBeta's native Ward-identity gauge check**, and the
  output records `native_rgbeta_checkinvariance_used=false`.
- The generated extra B/C(-1) `H Hdag S1dag S2dag` contraction was corrected
  from `eps(H,S2dag)*eps(Hdag,S1dag)` to
  `del(S2dag,H)*eps(Hdag,S1dag)` in the standard dual-doublet convention.
  A regression compares the two tensors under finite SU(2) rotations.
- The existing six B/C(-1) RGBeta couplings are **not modified**. Their
  independent exact tensor-normalisation checks remain outstanding.
- The runner retains the prior CLI fixes (`m1` for negative alpha; comma-separated
  model selection) and captures complete Wolfram logs.

## Local tests

    python -m pytest -q tests/python/test_general_t3_stage20.py tests/python/test_general_t3_stage20_quadratics.py

On the isolated packaging filesystem: 35 passed, 2 skipped (production-only
numerical-state tests). WolframScript is not installed in this environment;
no RGBeta beta-export success is claimed.

## Local Windows command

    powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage20.ps1 -ExportBetas -CheckGauge -Models B_-1,C_-1

The JSON output metadata is candidate-only. Do **not** use these couplings as
certified Matchete parameters or as complete physical neutrino-mass predictions.
