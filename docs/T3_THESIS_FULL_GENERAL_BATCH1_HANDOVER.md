# Batch 1 handover — A(-4), D(-2) full-inventory renormalisable running

**Scope:** UV one-loop renormalisable running, grouped F projection and EFT1 one-loop renormalisable running; not one-loop Wilson matching or the final neutrino mass evaluator. This batch is a numerical predecessor to Batch 2 and deliberately **does not** use the analytic hard-triangle/conditional restricted pilot matching.

## Files to install (repository-relative)

- `scripts/RunT3ThesisGeneralBatch1.py` — thesis-boundary runner, native RGBeta source inventory and F/EFT1 threshold handoffs.
- `Numerical/parallel_t3/background/general_batch1_api.py` — checksummed, read-only numerical interface to Batch 2; complex scalar/matrix values and interpolated trajectories.
- `tests/test_t3_thesis_general_batch1.py` — interface/contract smoke tests.
- `docs/T3_THESIS_FULL_GENERAL_BATCH1_HANDOVER.md` — this note.

No existing workstream, Python solver, Wolfram script, prior output or source config is replaced. The canonical general-interaction symbolic cache is the repository's existing `RunT3GeneralSymbolicStudy.py` pipeline.

## Installation and initial run (PowerShell)

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
Set-Location $repo

# The docs ZIP supplied earlier might not have been installed on this machine.
Test-Path .\docs\T3_A4_D2_GENERAL_PIPELINE_FEASIBILITY_2026_11_10.md
Test-Path .\docs\T3_A4_D2_FULL_GENERAL_IMPLEMENTATION_BATCHES.md
Test-Path .\configs\thesis_numerics\t3_thesis_A4_Dm2_pilot.json

# Replace $archive with the actual ZIP downloaded from the handover.
$archive = 'C:\Users\YOUR_USER\Downloads\T3_FullGeneral_Batch1_2026_11_10.zip'
Expand-Archive -LiteralPath $archive -DestinationPath $repo -Force

python -X utf8 -m unittest discover -s tests -p test_t3_thesis_general_batch1.py -v

# Requires wolframscript on PATH, Matchete and RGBeta installed.
# This creates/validates the full general symbolic UV/EFT1 exports once.
python -X utf8 -m scripts.RunT3ThesisGeneralBatch1 `
  --all `
  --symbolic-run-id thesis_full_general_sym_v1 `
  --run-id thesis_full_general_batch1_v1 `
  --generate-symbolic
```

If the **same** symbolic cache exists and is valid, omit `--generate-symbolic` for later Yukawa trials and use a **new** numerical `--run-id`. The underlying canonical symbolic tool uses file hashes to invalidate stale caches. Existing output runs are not overwritten. Use `--model A_-4` or `--model D_-2` instead of `--all` to diagnose one model.

If the symbolic exporter fails, first run its direct command to obtain its detailed Wolfram logs:

```powershell
python -X utf8 -m scripts.RunT3GeneralSymbolicStudy `
  --model A_-4 --run-id thesis_full_general_sym_v1
python -X utf8 -m scripts.RunT3GeneralSymbolicStudy `
  --model D_-2 --run-id thesis_full_general_sym_v1
```

The runner expects the original thesis `configs/thesis_numerics/t3_thesis_A4_Dm2_pilot.json`. It deliberately preserves the original BSM Yukawa matrices, `lambdaT3=0.1`, other scalar couplings `0.01`, and the thesis UV masses and gauge/Higgs inputs; SM flavour Yukawas use the existing Workstream-B initializer, explicitly recorded in the output. D(-2) `y2Conj` and extra holomorphic scalar mass have **zero UV boundary values** but their beta functions are still integrated; the five extra scalar quartics start at `0.01` and run.

## Numerical output structure

For each model `A_-4` and `D_-2`:

```text
output/parallel_t3/full_general_batch1/<run-id>/<model>/
    manifest.json
    trajectories.npz
    fermion_matching_input.json
    scalar_matching_input.json
```

A shared `output/parallel_t3/full_general_batch1/<run-id>/summary.json` records success or the explicit blocker.

- `manifest.json` records the **full native UV and EFT1 beta inventories**, every complex-valued endpoint, UV/F/S scales, scalar extra-slot names, source SHA256, input benchmark, numerical solver counts, optional explicit F-boundary shifts, and checksums of the three output artifacts.
- `trajectories.npz` has `uv__mu_gev`, `eft1__mu_gev` and `uv__<beta symbol>` / `eft1__<beta symbol>` for **every** native report beta. Scalar values are complex128 vectors over scale; matrix couplings are complex128 `(N,3,3)` arrays. The complete original UV Yukawa matrices and their run to F are present. No beta keys are silently discarded. Every scale sample has the full native parameter state.
- `fermion_matching_input.json` contains full UV state at `muF=1e5 GeV`, F masses and Yukawas (`MF`, `y1`, `y2`, and `y2Conj` in D(-2)), as well as the complete native EFT1 tree-projected boundary. It explicitly requests removal of **all F at the same scale**; it does **not** assert that the one-loop Wilson matching has happened.
- `scalar_matching_input.json` contains full scalar-only EFT1 state at `muS=1e3 GeV`, both scalar masses, all scalar quartics including D(-2)'s five extra native slots, and D(-2)'s complex `muT3GeneralS2Sq`. It also includes the existing **candidate** real-scalar quadratic Hessian in the spherical basis. The Matchete/RGBeta mass-basis map is explicitly not claimed to be complete.

## Batch 2 Python interface

```python
from Numerical.parallel_t3.background.general_batch1_api import load_batch1

run = load_batch1(
    'output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1/D_-2'
)

f = run.fermion_matching_input()  # JSON, all fermions still active above F
m = run.at_fermion_threshold()     # native numpy values at F
print(m['MF'].shape, m['y1'].shape, m['y2'].shape, m['y2Conj'].shape)

# SM + BOTH scalar multiplets, F removed; all general scalar terms retained.
s = run.scalar_matching_input()
eft1_1e4 = run.at('eft1', 1e4)       # log(mu)-interpolated background
end = run.at_scalar_threshold()    # exact saved numerical endpoint at S
print(end['mS1Sq'], end['mS2Sq'], end['muT3GeneralS2Sq'])
```

**Essential for Batch 2:** the symbols and matrices in these objects are *native RGBeta basis coordinates*, not normalized Matchete Wilson operators. The one-loop matcher must use `Lagrangian/RunMatching.wl::RunSequentialT3Matching` and derive or consume the explicit Matchete↔RGBeta coupling and scalar mass-basis transformations. Match **all F** at F and run generated EFT1 Wilson coefficients alongside this background before matching **both S** at S. The complete one-loop F matching result must include its UV/EFT1 pole and finite subtraction prescription. It must not be replaced with the analytic hard triangle.

### Renormalisable one-loop threshold corrections (future Batch 2)

The default F boundary is a transparent *tree-level projection for renormalisable parameters*, not a claim that one-loop finite parameter shifts vanish. If Batch 2 derives local renormalisable shifts, save an **additive** JSON file with exact native names and re-run EFT1:

```json
{
  "format": "t3-batch1-eft1-shifts-v1",
  "model": "D_-2",
  "mu_F_GeV": 100000,
  "deltas": { "lambdaH": 0.001 }
}
```

The displayed `0.001` is **format illustration only**, not a calculated T3 threshold correction. Apply a real calculated file using `--eft1-boundary-shifts path/to/shifts.json` with a fresh `--run-id`. The runner refuses unknown couplings and re-integrates EFT1, while keeping the exact original UV/F state. The manifest and trajectory loader record and validate every declared discontinuity. This gives Batch 2 a real re-integration interface when the matching produces renormalisable boundary corrections.

## Status / accuracy

This batch reuses the repository's existing `UVRunner.py`, `IntermediateScalarRunner.py`, `GeneralRunningBenchmark.initial_state`, beta-expression evaluators, scalar quadratic Hessian and tree-level F projection. It **does not** redefine the beta functions. No one-loop Wilson coefficients are fabricated. The source RGBeta export itself is candidate/native and the physical matching dictionary is still an open implementation step.

**Local smoke tests:** imports, complex matrix encode/decode, native D(-2) extra-key endpoint storage, checksum protection, exact boundary-shift handling, interpolation, and rejection of a missing `y2Conj` operator inventory. These tests do **not** demonstrate a fresh Wolfram export or source-RGE numerical solve; the latter must be executed on the Windows machine with the actual repository and Wolfram installations.
