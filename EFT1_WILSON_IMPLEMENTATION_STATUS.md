# EFT1 Wilson ODE — scoped implementation status

GitHub source reviewed (default branch, 2026-11-09):

- `RGE/running/intermediate/GeneralEFT1C5Mixing.py` — seven-term `calculate_master_beta` for master-basis `psi^2 phi^2`, source projectors, and fixed-order `solve_boundary`.
- `Numerical/running/GeneralT3SequentialLeadingLog.py` — conditional scalar-pair `LLSS -> C5` leading log and finite bubble.
- `Numerical/running/GeneralT3DynamicFermionWilson.py` — dynamic nonmixed F-tree tensors.
- `Numerical/running/GeneralT3MixedExactFermionWilson.py` — native mixed B(-1)/C(-1) F-tree tensor.
- `Numerical/orchestration/GeneralRunningBenchmark.py` — UV and EFT1 renormalisable endpoint snapshots only, though internal EFT1 runner retains a trajectory.
- `Numerical/running/IntermediateScalarRunner.py` — uses `solve_ivp`, `t=ln(mu)` and stores the intermediate coupling trajectory.
- `Numerical/joint_scalar_matching/large_batch.py` — `C5` fixed-order log and joint finite bubble once, explicitly without `LLSS` self-running.
- `scripts/RunT3All16DynamicSequential.py`, `scripts/RunT3All16DynamicSequential16.py`, `scripts/RunT3SequentialThesisRefresh.py` — current sequential orchestration and donor fits (files present on GitHub at inspection time).

## Generated paths

- `Numerical/running/GeneralEFT1WilsonODE.py` — **new**, independent solver for the already implemented seven-term master `psi^2 phi^2` anomalous dimension, now with dynamically evolving Wilson tensor and interpolated renormalisable coupling trajectory.
- `scripts/TestEFT1WilsonODE.py` — **new**, small synthetic smoke test with two retained scalar blocks.

## What is implemented

`solve_coupled_wilson(boundary, snapshots, ...)` integrates the entire symmetric `C_ijab` coefficient array using the repository's existing Eq. (4.85) evaluator. Every coupling snapshot is a complete canonical-Weyl payload accepted by the existing `load_boundary`, plus `mu_gev`; snapshots must span `mu_S` and `mu_F`. Gauge representation generators must remain identical. All field-index tensor components are evolved, including induced `LLHH`. The latter is reported by `project_c5`. The new solver does **not** independently add the old leading logarithm.

## What is NOT implemented

This is *not yet connected to the real T3 sequential runner*. The available Matchete F-tree tensor has shape `(3,3,2,2,n,n)` and a **formal ordered-spinor convention**, while the master RGE requires fully normalized Weyl and real-scalar component tensors `(nf,nf,ns,ns)`. The source itself labels the signed `GammaCC -> Weyl` conversion unverified. The present `scalar_tensors` returns only `lambda_HHSS`, not the **entire** `lambda_abcd` required for the master self-running. The master evaluator further requires a complete set of gauge generator matrices and light-field Yukawa tensors. The existing orchestration JSON also discards EFT1 trajectory samples. None of these conversions can be replaced by assumed identities or zeros without changing the claimed physics.

In addition to Eq. (4.85), the cited general EFT paper provides the `psi^2 phi^2` mixing from additional operator classes in Eqs. (4.84), (4.86)–(4.91) and beyond. The complete source-operator inventory has not been demonstrated for all T3 models. In particular, these cross-class mixings must be considered if the corresponding operators have nonzero leading matching coefficients. This new solver integrates **only** the Eq. (4.85) block, including its self-running and internal mixing.

Therefore the new files are deliberately **additive**. They do not alter the historical matching, optimisations, or thesis plots, and do not mislabel a conditional result as physical.

## PowerShell

Extract the ZIP into the repository root, keeping directories:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m scripts.TestEFT1WilsonODE
```

To integrate a **previously certified** master-basis boundary and matching coupling trajectory:

```powershell
python -X utf8 -m Numerical.running.GeneralEFT1WilsonODE `
  --boundary input\eft1_canonical_master_boundary.json `
  --coupling-trajectory input\eft1_canonical_coupling_trajectory.json `
  --output output\general_sequential_pilot\eft1_master_ode_v1\D_m2.json
```

For synthetic examples only, append `--allow-synthetic`. These canonical inputs **do not yet exist** for the physical D(-2)/A(-4) pipeline; the command is an interface example, **not** a validated benchmark run instruction.

## Testing performed in this environment

- Python syntax compilation: successful.
- Synthetic smoke test of solver logic with a minimal local stand-in for the existing `calculate_master_beta` **scalar-pair-only, vanishing-Yukawa/gauge** case: successful, `nfev = 629`.
- Synthetic result: `C5[0,0] numerical = -1.163109663522e-7 GeV^-1`; fixed-order leading log `-1.166503235297e-7 GeV^-1`.
- **Not performed:** import against the checked-out full repository, full seven-term physical benchmark, D(-2)/A(-4), any 16-model run. Do not interpret this smoke test as physics verification.

## Next required input

For a real local pipeline integration, the most useful input is a **ZIP of the current local repository** or at least the local EFT1 Matchete/RGBeta operator exports, scalar quartic dictionary, and exact sequential entrypoint files. The source must be checked for the signed Weyl matching convention before a physical RGE can be switched on.
