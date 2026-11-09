# EFT1 Stage B: A(−4) same-benchmark numerical transport

**Scope:** Additive runner for the *existing* A(−4) benchmark. No GitHub modifications and no replacement of the validated fixed-order calculation. This does **not** claim physical Wilson normalization.

## What is newly implemented

- Load the repository's real A(−4) sequential `dynamic_matching.npz`, which already contains the original 3-flavour matched `C_LLSS_muF` and per-operator sources.
- Verify the source tensor's shape `(3,3,2,2,8,8)`, complex tensor exchange symmetries, source decomposition, scale metadata and **SHA-256 of the actual Matchete seed bytes plus tensor bytes**, against the saved sequential `result.json`.
- Read the same run's actual `renormalisable_parameters_at_F` rather than a different UV fit or the unrelated `output/a_minus4_eft1_batch2` tensor snapshot.
- Re-integrate the **already archived** A(−4) 17-coupling RGBeta EFT1 beta functions for 7 points from `M_F` to `M_S` with both scalars active. Confirm the reconstructed endpoint matches the saved benchmark `renormalisable_parameters_at_S`. If not, fail rather than joining incompatible files.
- Build the 45-Weyl / 12-real-scalar running Yukawa and quartic tensors and electroweak/strong generators from the checked `Aminus4EFT1Couplings` module. Run the previously installed Stage A `T3EFT1WilsonODETransport` in explicitly **formal identity-map** mode.

**Existing input paths** (checked on public `main` branch on 2026-11-09):

1. `output/general_sequential_pilot/all16_dynamic_scalar_threshold/batch16_v1/A_m4/result.json`
2. `output/general_sequential_pilot/all16_dynamic_scalar_threshold/batch16_v1/A_m4/dynamic_matching.npz`
3. `output/general_observable/fermion_first_all16_v1/A_m4/data/eft1_after_F_wilson_seed.json`
4. `output/studies/general_scalar_extended/complete_scalar_symbolic/full16_general_v1/A_m4/eft1_rgbeta.json`

The saved benchmark's hierarchy is `mu_UV=1e7 GeV`, `mu_F=1e5 GeV`, `mu_S=1e3 GeV`, `mu_low=1e2 GeV`. Values are **read and validated from the saved inputs**, not hard-coded into the runner.

### Important source discrepancy avoided

The separate `output/a_minus4_eft1_batch2/A_m4_eft1_tensors.npz` starts at e.g. `g2(M_F)=0.65404`, while the saved `batch16_v1` dynamic sequential benchmark has `g2(M_F)=0.59802`. They are **not interchangeable**. This runner never borrows the first trajectory. It instead reconstructs its EFT1 tensors from the saved sequential benchmark's own `M_F` endpoint and checks its `M_S` endpoint.

## Files to install

- `scripts/RunAminus4SavedEFT1Pilot.py`
- `tests/test_a_minus4_saved_eft1_pilot.py`
- `EFT1_STAGE_B_SAVED_BENCHMARK.md`

These three files are additive. Stage A's two modules must already be installed at:

- `Numerical/running/T3EFT1MasterSourceBridge.py`
- `Numerical/running/T3EFT1WilsonODETransport.py`

## PowerShell commands

From `C:\thesistemp\T3neutrinomassmodelpipeline` with the existing `.venv` activated:

```powershell
$zip = Join-Path $HOME 'Downloads\T3_EFT1_StageB_Saved_Benchmark.zip'
Expand-Archive -LiteralPath $zip -DestinationPath .

python -X utf8 -m unittest discover -s tests -p 'test_a_minus4_saved_eft1_pilot.py' -v

# First verify that the saved fermion matching source and rebuilt RGBeta EFT1
# running trajectory are actually from the SAME saved benchmark.
python -X utf8 -m scripts.RunAminus4SavedEFT1Pilot --prepare-only

# Once the above passes, use a different output directory for the full ODE:
python -X utf8 -m scripts.RunAminus4SavedEFT1Pilot `
  --output 'output\eft1_stageB_saved_a_m4\batch16_v1_full_ode'
```

Each invocation refuses to overwrite an output directory. The `--prepare-only` outputs go under `output/eft1_stageB_saved_a_m4/batch16_v1/` by default.

For a quick single-stage run (without `--prepare-only`) in a new directory:

```powershell
python -X utf8 -m scripts.RunAminus4SavedEFT1Pilot `
  --output 'output\eft1_stageB_saved_a_m4\initial_ode'
```

**Generated output locations:**

- `output/eft1_stageB_saved_a_m4/<run>/A_m4_formal_F_source.npz` and `.json` (verified original matched source; spinor normalization unresolved)
- `output/eft1_stageB_saved_a_m4/<run>/A_m4_eft1_tensors.npz` and `A_m4_eft1_tensor_basis.json` (matched native running input)
- `output/eft1_stageB_saved_a_m4/<run>/input_validation.json` (matching hash and scalar running endpoint check)
- `output/eft1_stageB_saved_a_m4/<run>/pilot_summary.json` (formal ODE C5 projection, Ward and symmetry residuals)
- `output/eft1_stageB_saved_a_m4/<run>/formal_master_C_at_muS.npz` (full formal `C_ijab` at `M_S`)

## Scientific limits and required next physics work

1. The saved source and its 3×3 flavour matrices are real numerical inputs, but the overall `GammaCC` → master-Weyl sign/factor remains **unproven**. The pilot uses factor `1` strictly for illustration.
2. The native Matchete CG → RGBeta scalar basis and quartic coupling normalization remain **uncertified**. The pilot uses an explicitly enabled *formal* identity coordinate map. This includes a previously documented A(−4) portal-normalization mismatch in the RGBeta reference-basis quartic audit (`lambdaH2Adj` coefficient `3` versus `1`).
3. The seven-term master equation evolves the `ψ²φ²` class only. Closure with any other applicable dimension-five classes must be established separately before a physical result.
4. Do **not** add the old fixed-order `C5_direct_log` or scalar finite bubble to the evolved result. At the eventual scalar threshold, subtract the shared infrared logarithm **once** and add the separately matched, nonduplicated finite contribution. Until that point neither neutrino masses nor fit values from this run are physical.
5. This script has been syntax-checked and tested against synthetic saved-file fixtures; the user's local `batch16_v1` files and full ODE have **not** been run in this environment.
