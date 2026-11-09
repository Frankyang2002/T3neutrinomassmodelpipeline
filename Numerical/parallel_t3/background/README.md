# Workstream B — all-16 renormalisable background and formal all-F boundary

**Package revision:** B1 (9 November 2026). **Read-only GitHub main inspected:**
`a944e0030b8ef24b23342d99e0d8f33cef3cffd5`.

This addition **does not edit** `Lagrangian/`, `RGE/`, the symbolic exports, Workstream A/C
or any GitHub file. The user's local checkout and the Batch 1 native outputs are
required for actual physics integration. This ZIP includes the executable B code and
tests, **not fake precomputed trajectories**. Our isolated container has neither
that checkout nor Wolfram/Matchete. Run the PowerShell commands below on the local machine.

## Exact PowerShell instructions

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
# Extract the supplied ZIP at the project root. The archive already contains
# Numerical/parallel_t3/background/, scripts/, and tests/ path prefixes.
Expand-Archive -LiteralPath 'C:\PATH\TO\T3_Workstream_B_all16.zip' -DestinationPath . -Force

# Existing .venv is expected to contain numpy, scipy and pytest.
.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_parallel_b_contracts.py
$env:T3_B_RUN_NATIVE_TESTS='1'
.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_parallel_b_native.py
Remove-Item Env:T3_B_RUN_NATIVE_TESTS

# First test A(-4) against the *corrected* Matchete seed.
.\.venv\Scripts\python.exe -X utf8 scripts\RunT3ParallelBBackground.py `
  --model A_-4 --run-id B_A4_trial_01

# Run the complete 16-model batch (reuse existing Batch 1 exports).
.\.venv\Scripts\python.exe -X utf8 scripts\RunT3ParallelBBackground.py `
  --all --run-id B_all16_qualitative_v1 `
  --mu-uv 1000000 --mu-f 100000 --mu-s 1000 `
  --uv-samples 17 --eft1-samples 33

# Display the per-model dispositions; blocked matching does not erase RGE data.
$s = Get-Content 'output\parallel_t3\background\B_all16_qualitative_v1\summary.json' -Raw | ConvertFrom-Json
$s.model_results.PSObject.Properties | ForEach-Object {
    [pscustomobject]@{Model=$_.Name; Background=$_.Value.background;
      FTree=$_.Value.fermion_boundary; HardC5=$_.Value.hard_C5_status;
      Blocker=$_.Value.blocker}
} | Format-Table -AutoSize
```

If a model was already attempted, use a **new run-id**. `--force` permits
recomputing that B-owned run-id; it clears only that model's B output files,
not the source or Batch 1 outputs. Changing `--mu-f` away from 100000 GeV
currently blocks the archived B_-1/C_-1 mixed F-tree seed, which hardcodes
its own matching scale; ordinary 14 models are not hardcoded.

### Inputs: explicit contract

| Input | Required provenance / semantics |
|---|---|
| `output/batch1_general_contract/batch1_native_v1/models/<tag>.json` | Source SHA256 per file, UV/EFT1 cached-output hashes, and F inventory |
| `.../rgbeta/<tag>/uv.json` | `SuccessCandidateMassBeta`, full `report_betas`, UV metadata and `beta_symbols` |
| `.../rgbeta/<tag>/eft1.json` | Same, `GeneralScalarStage=EFT1`, no `y1/y2/y1Conj/y2Conj/MF` |
| `.../matchete/<tag>/data/eft1_after_F_wilson_seed.json` | For 13 non-A4 nonmixed points; actual source factorisation, CGRegistry, `IntegratedFields=[F]` |
| `output/general_observable/fermion_first_hermitian_portal_v1/A_m4/data/eft1_after_F_wilson_seed.json` | Corrected A(-4) source, SHA and both portal checks must agree with Batch 1 |
| `output/general_observable/ordered_wilson_tensor_v4/{B_m1,C_m1}/` | Dedicated previously derived mixed X/Y native source used only for B_-1/C_-1; source residuals and CGs checked in native adapter |

Root paths above are relative to
`C:\thesistemp\T3neutrinomassmodelpipeline`. No symbolic work is repeated.
Checksums are **byte SHA256**, so if Batch 1 artifacts were converted/reformatted,
regenerate their provenance consistently rather than silently accepting an old beta.
For the core `Lagrangian` and `RGE` files that produced Batch 1, a mismatch
blocks reuse. Changes to A-owned Wilson-only sources are recorded but not a
B beta-cache blocker. Public Git SHA is **not** used to assert local equality.

`GeneralRunningBenchmark.initial_state` is reused. All beta symbols must be
mapped to a known numerical type using the existing
`Numerical/core/GeneralCouplingLayout.spec_for`; no unknown name gets zero beta.
The specific benchmark overrides initialise distinct complex, nondiagonal
3x3 Yukawa matrices, diagonal initial heavy masses at 0.91, 1, 1.09 times MF,
nonzero invariant quartics and (where present) complex mixed/self quadratic
couplings. Remaining baseline parameters are supplied by the existing general
benchmark initializer. **The complete resolved initial UV state** (not just
those overrides) is saved in each manifest. This benchmark is qualitative,
not a neutrino/CKM fit or certified vacuum-stable parameter point.

### Actual numerical evolution

For each model, both segments integrate every beta in the source's native
complete inventory with the existing DOP853 runners. RGBeta convention:

\[
16\pi^2\frac{dX}{d\ln\mu}=\beta^{(1)}_X,
\qquad \mu_{UV}>M_F>M_S.
\]

- **UV**: all 3 gauge couplings, full 3-generation complex SM Yukawas,
  y1/y2/(conditional y1Conj/y2Conj), 3x3 complex heavy MF, all invariant
  quartics and allowed scalar quadratic/mixing structures.
- **F threshold at one MF**: all three heavy fermions are removed together,
  with the common SM+S1+S2 renormalisable coordinates copied from the actual
  evolved UV endpoint. This is **tree-level renormalisable projection**;
  finite one-loop renormalisable threshold shifts are **not calculated**.
- **EFT1**: the complete native scalar-only inventory, **both S1 and S2**
  retained; no y1/y2/MF. Save the numerical scalar couplings at MS, but do
  not integrate out either scalar in B.

The UV and EFT1 beta inventories are exact complements after removing the five
possible F coordinates. Every non-F value is numerically checked continuous
at MF. Nonfinite integrals and couplings above 4pi are rejected (rough veto,
not a loop-level theoretical perturbativity theorem).

### Output schemas, per model

Directory: `output/parallel_t3/background/<run_id>/<tag>/`.

| File | Content |
|---|---|
| `background.npz` | `mu_uv_gev[n_uv]`, `mu_eft1_gev[n_eft]` descending, `t_uv`, `t_eft1`; **every** `uv__<beta_name>` and `eft1__<beta_name>`. Each coupling trajectory is `complex128` of shape `[samples]` or `[samples,3,3]`; scales are `float64`. No pickles. |
| `background_manifest.json` | exact model and scales, both ordered beta inventories, all resolved UV starting values, array shapes, numerical solver diagnostics, source and export SHA256s, local git status, SHA256 of NPZ, full native scalar coordinates at MS, certification flags |
| `fermion_boundary_formal.npz` | **only if available**: `C_LLSS_formal[3,3,2,2,n,n]`, per-source `C_LLSS_by_source_formal[k,3,3,2,2,n,n]` in source Matchete commuting-spinor coordinates, `source_indices[k]`, `heavy_mass_eigenvalues_gev[3]`, with `n=2(dS1+dS2)`. Optional `C5_hard_external[3,3]` if and only if externally calculated. |
| `fermion_boundary_manifest.json` | source and signature SHA, exact named axes and heavy/flavour basis, full hard-C5 status, one-loop incompleteness, special-mass blocker/eigenbasis status |
| `blocker.json` | on failed renormalisable stage, precise exception and stack tail; `summary.json` always includes all models attempted |

The **formal** LLSS coefficients are in a Matchete source polynomial: no
GammaCC↔Weyl factor is implied, and their absolute physical coefficient is
**not** certified. Same-run UV Yukawas/MF are factored into **every** source
operator via the existing `flavour_matrix`, SVD/Takagi, CG and orbit routines,
not a guessed `y1/MF@y2.T` substitute. Its full complex 3x3 flavour indices
are preserved; the scalar indices follow the native real scalar field basis.

For B_-1/C_-1 the existing `GeneralT3MixedExactFermionWilson` native X/Y
interface is used with the **MF** scalar quadratic state, **not** with a
restricted scalar potential. This is a candidate because the mixed X/Y source
and RGBeta mass basis are not independently identified. If native sources,
residual audits or the scale disagree, only the F boundary is marked blocked;
its successful renormalisable trajectory is retained.

Six models `A_-2,A_+0,D_-2,D_+0,E_-2,E_+0` have neutral self-bilinear
complex-phase conventions **not certified** by the existing native real-split
mass code. B saves their complex coupling and does **not** force it to a real
number. Workstream C must resolve the phase/eigenbasis identification. This
may block subsequent physical scalar matching even if B's formal LLSS tensor
is available.

### Read-only cross-workstream API

```python
from pathlib import Path
from Numerical.parallel_t3.background.api import Background, load_formal_boundary
folder = Path('output/parallel_t3/background/B_all16_qualitative_v1/A_m4')
bg = Background.load(folder)  # Verifies NPZ SHA, full beta inventory and MF continuity
couplings_mf_to_ms = bg.native_eft1_at(1.0e4) # 3x3 complex SM Yukawas included
scalar_coords_ms = bg.native_scalar_ms()       # native invariant couplings / quadratic coordinates
formal, formal_meta = load_formal_boundary(folder) # exact axis labels, source hash
```

These are **native RGBeta** coupling coordinates and **formal Matchete** LLSS
coefficients, respectively. Workstream A must provide the explicit signed
operator map and full coupled dimension-five ADM. Workstream C must construct
mass/eigenbasis and quartic component tensors at MS, and match both scalars
**simultaneously**. B supplies neither an artificial C5(MF)=0 nor a second
fixed-order LLSS→C5 leading-log calculation.

### Optional independently calculated nonzero one-loop hard C5

Use `--hard-matching-root <relative_dir>` only for an independently computed
same-run fermion-threshold matching result. Each model's
`<relative_dir>/<tag>/hard_C5.npz` must contain `C5_hard[3,3]` as `complex128`;
its JSON companion `hard_C5_manifest.json` must have:

```json
{
  "model": "A_-4",
  "mu_F_gev": 100000,
  "same_run_uv_state_sha256": "<copied from background_manifest.UV_endpoint_sha256>",
  "method": "independent_one_loop_fermion_threshold_calculation",
  "basis": "<explicit operator, flavour, phase, normalization convention>",
  "hard_c5_file_sha256": "<SHA256 of hard_C5.npz>"
}
```

This is an interface contract, **not a claim that hard matching exists**.
A missing or unverified hard term is recorded as unavailable/blocked; no
physical zero is supplied. Computing one hard C5 alone still does not certify
all other one-loop hard dimension-five terms, matching completeness, or basis
conversion.

### All-16 work readiness (not fabricated run outcomes)

| Model | Candidate beta pair (Batch 1) | Formal F-tree approach | Special limitation |
|---|---|---|---|
| A_-4 | UV+EFT1 | Corrected source/CG | Anti-Hermitian portal must remain fixed |
| A_-2, A_+0 | UV+EFT1 | 14-model nonmixed native | Neutral self-bilinear complex phase uncertified |
| B_-3, B_+1 | UV+EFT1 | 14-model nonmixed native | Matchete/RGBeta dictionary unresolved |
| B_-1 | UV+EFT1 | Dedicated native X/Y | Mixed mass eigenbasis mapping candidate |
| C_-3, C_+1 | UV+EFT1 | 14-model nonmixed native | Matchete/RGBeta dictionary unresolved |
| C_-1 | UV+EFT1 | Dedicated native X/Y | Mixed mass eigenbasis mapping candidate |
| D_-2, D_+0 | UV+EFT1 | 14-model nonmixed native | Neutral self-bilinear complex phase uncertified |
| D_+2 | UV+EFT1 | 14-model nonmixed native | Matchete/RGBeta dictionary unresolved |
| E_-4, E_+2 | UV+EFT1 | 14-model nonmixed native | Matchete/RGBeta dictionary unresolved |
| E_-2, E_+0 | UV+EFT1 | 14-model nonmixed native | Neutral self-bilinear complex phase uncertified |

Actual numerical successes/failures are recorded **per model only after**
executing this runner on the local repository. No workstream-A ADM or C scalar
matching is claimed. Tests check basic source/beta/hash contracts and a mock
array roundtrip; they do not certify the full physical normalization.
