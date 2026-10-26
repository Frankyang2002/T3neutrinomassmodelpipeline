# Stage 21: 16-model candidate RGBeta export validation

Stage 20 produced 32 Wolfram RGBeta JSON exports: UV and scalar-only EFT1 for each of the 16 ordinary T3 models. Stage 21 is a **read-only compatibility and coverage audit** of those outputs.

## What is checked

- All 16 registry assignments, 50 model-specific extra scalar-quartic slots: 38 new RGBeta vertices and 12 structurally reused B/C(-1) slots.
- All 32 exported JSON payloads are successful, have matching representation and stage metadata, and share the model-specific registry beta keys.
- Both newly generated and previously existing overlapping quartics are included in RGBeta's **full quartic projection** (`BetaTerm[Quartic,1]`).
- Every serialized beta is a nonempty string and no unresolved RGBeta beta-projection heads remain.
- UV/EFT1 beta key sets remain consistent apart from Yukawa/heavy-fermion fields absent after integrating out F.
- The report leaves exact Matchete normalisation, quadratic mixing, and physical C5 matching explicitly **unverified**.

`--numeric-probe` separately uses synthetic complex Yukawa matrices (3x3) to exercise the production Python InputForm expression evaluator. A failure is a diagnostic for parser or synthetic input shape problems; passing is **not proof** of a representation-correct numerical RGE trajectory. Examine each failure in `output/general_t3_stage21/stage21_summary.json`; prefer actual baseline state tests next.

## Install

Unzip the Stage21 overlay at the project root, preserving paths.

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage21_ExportAudit.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_stage21.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage21.ps1
```

Optional diagnostic:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage21.ps1 -NumericProbe
```

Required input files (already created by Stage 20): `RGE/general/generated/general_rgbeta_quartic_registry.json` and `output/general_t3_stage20/{MODEL}_{uv,eft1}_rgbeta.json`.

Output: `output/general_t3_stage21/stage21_summary.json`.

## Physical status

The verified result, if all checks pass, is *candidate quartic-RGE export consistency*. It is **not** certified general scalar potential matching or a complete neutrino mass prediction. Next steps: confirm Matchete↔RGBeta invariant transforms and signs, run model-specific state-vector ODE tests, implement eight scalar bilinear mass mixings, and recompute masses/thresholds in the correct basis.
