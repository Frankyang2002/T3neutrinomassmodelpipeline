# Worker B next: running source bridge (A_-4 / D_-2)

## Status
Standalone bridge, no production changes. Uses saved Batch3B cached affine full-real quartic response, replaces its HHSS orbit with caller-supplied source HHSS, and constructs *unbroken* scalar Hessian from running EFT1 masses. Requires the corresponding cached source files; do not mistake unit tests for full local benchmark validation.

## Files and destinations
- `Numerical/running/T3GeneralScalarBackground.py`: previous complete module, included to make ZIP standalone.
- `Numerical/running/T3ScalarRunningSourceBridge.py`: new running-state and candidate quartic adapter; D_-2 native complex bilinear.
- `scripts/CheckT3WorkerBScalarSource.py`: read-only endpoint runner, requires same-scale source HHSS exported by existing project code as NPZ key `hhss`.
- `tests/test_t3_scalar_running_source_bridge.py`: synthetic mass/affine/bridge checks.
- `docs/WORKER_B_NEXT_README.md`: this file.

## Source conventions
- Public reference files inspected: `Numerical/running/T3CoupledQuarticContext.py`, `Numerical/joint_scalar_matching/special_mass_v6.py`, `Numerical/joint_scalar_matching/numerical_benchmark_v5.py`.
- Source `special_mass_v6.py` explicitly defines `L_int=+mu (S_dag)^2+h.c.` for D_-2 and a native real quadratic block `[[m2-2 Re mu,-2 Im mu],[-2 Im mu,m2+2 Re mu]]`.
- A_-4: diagonal 2+6 real BSM block; D_-2: diagonal 6 S1 plus complex split 2 S2 block.
- The Higgs mass block is set to zero because this interface exposes only the heavy scalar unbroken Hessian. Do not use it for physical Higgs masses.
- Batch3B's non-HHSS basis remains *candidate*; only externally injected source HHSS is source-based. No claim of complete physical potential or final matching.
- No sourced all-16 RG quartic realization or B(-1)/C(-1) mixings are claimed.

## PowerShell
```
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive "$HOME\Downloads\T3_Worker_B_Running_Source_Bridge.zip" -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p 'test_t3_scalar_running_source_bridge.py' -v
```

If the existing archived Batch3B outputs and independently evaluated source HHSS arrays are present, use:
```
python -X utf8 -m scripts.CheckT3WorkerBScalarSource --model A_-4 `
  --response-root 'output\audits\batch3b_quartic_fixed_20261111_234059' `
  --hhss-npz 'PATH_TO_SAME_SCALE_A4_SOURCE_HHSS.npz' `
  --output 'output\audits\worker_b_source_a4_20261112_v1.json'
```
Repeat for D_-2 with its own same-scale source HHSS. The runner intentionally fails if the HHSS file is absent or the Batch3B provenance hashes mismatch. Do not use the illustrative placeholder HHSS path literally.

## Remaining blockers
1. Local-source all-16 quartic inventory mapping and native/RGBeta absolute normalizations.
2. Same-scale source HHSS supplied for both benchmark models (not present in the uploaded ZIP).
3. B/C dedicated mixing source integration, all additional model-specific bilinears.
4. Independent comparison against saved local physical/candidate scalar benchmarks.
