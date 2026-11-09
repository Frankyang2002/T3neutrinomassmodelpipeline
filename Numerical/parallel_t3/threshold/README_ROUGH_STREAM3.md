# T3 Parallel Stream 3 — rough fermion-first all-16 (2026-11-10)

**Source inspected read-only:** `Frankyang2002/T3neutrinomassmodelpipeline` main
at `25eb5e13d09dd2f138176d4e36c4fc65edfa4bb6` (GitHub). Uses the existing
`threshold/matching.py`, `threshold/pipeline.py`, `threshold/preflight.py`,
`threshold/models.py`, and `Numerical/running/SMWeinbergEvolution.py`. No GitHub
writes or symbolic regeneration. The user's *uncommitted* local files could not
be inspected; check for local conflicts before extracting.

## EFT and numerical scope

`UV → all F at one MF → EFT1 with both S1/S2 and the complete Stream-2 Wilson
endpoint → all scalar eigenstates at one MS=1000 GeV → SMEFT → SM+Weinberg
one-loop running to 173 GeV → neutrino mass/mixing diagnostics.`

The matching calculation is **exactly the pre-existing C scalar-threshold
solver**, evaluated with `mode='rge_improved'`. The pre-existing C solver
computes `C5_SMEFT(MS) = C5_EFT1(MS) + ΔC5_scalar_bubble(MS) +
C5_independent_finite(MS)` and includes *all ordered eigenstate pairs* (or all
unordered pairs with multiplicity) with the finite mass-dependent `B0bar`.
It does NOT add Stream 2's LLSS→LLHH leading log again. When the independent
hard term is missing, the existing solver sets that term to zero; this is
labelled an explicit **rough benchmark approximation**, not evidence the term
vanishes. Other independent finite scalar/wavefunction/EOM terms are uncomputed.

The existing C runner diagonalises the full real symmetric scalar Hessian. A
nonzero complex B/C mixed-doublet bilinear or neutral self-bilinear appears as
real R–I off-diagonal elements in the supplied Hessian, not as an imaginary
entry to the real Hessian. The bridge *does not discard* those elements. All
heavy scalar eigenstates are integrated **at the same MS**, even if their
masses differ. The optional `models.from_parameters` source-based constructor
is only used if the B owner explicitly declares
`mass_hessian_constructor='explicit_t3_native_quadratic_v1'` and a canonical
scalar frame; otherwise a complete B Hessian is mandatory.

The SM ODE is the original `SMWeinbergEvolution.evolve_sm_weinberg`.
The Majorana mass matrix uses `mν=-(246.22 GeV)^2 C5/2` and sorted Takagi
singular values. PMNS mixing uses the **left** singular vectors of the
charged-lepton `ye(MS)`; those vectors are invariant under the existing
one-loop SM-only `beta_ye = a ye + 3/2 ye ye† ye`. Angles are omitted if
Takagi or charged-lepton singular vectors are numerically degenerate. Sorted
column order is not automatically normal/inverted phenomenological ordering.

## Producer input directories (exact handoff)

- **Stream 1**: `output/parallel_t3/rough_inputs/rough_fermion_first_v1/c/<slug>/background.json` and `background.npz`.
- **Stream 2**: `output/parallel_t3/rough_wilson/rough_fermion_first_v1/<slug>/wilson.json` and `wilson_at_MS.npz`.
- `<slug>` is `A_m4`, `A_m2`, `A_p0`, `B_m3`, `B_m1`, `B_p1`, ..., `E_p2`.

Supported input semantics:

1. The complete published **Workstream C** `t3-parallel-c-wilson-input-v1`
   and `t3-parallel-c-background-input-v1` JSON/NPZ pairs are accepted and
   verified, not renamed blindly.
2. The existing **Workstream A** `t3-parallel-a-wilson-output-v1` full `C_full`
   tensor is converted using explicit 3-family Weyl/scalar label lists and
   the declared no-factorial Weinberg normalization. A legacy `manifest.json`
   is discoverable in place of `wilson.json`, but must still contain current
   explicit matching state IDs and the synthetic/genuine flag.
3. The existing **Workstream B** `t3-parallel-B-background-v1` can supply actual
   `eft1__gY`, `eft1__g2`, `eft1__g3`, `eft1__lambdaH`, `eft1__ye`,
   `eft1__yu`, `eft1__yd` sampled at the saved MS. Those are converted to
   the *true final time slice*, checking scale ordering, complex matrix shapes,
   and numeric finiteness. Its native `eft1__lambda...` coordinates do **not**
   suffice to reconstruct a full physical HHSS quartic. It therefore requires
   an additional source-normalised, *fully symmetric* `quartic_real_full_ms`
   tensor and a full `mass_squared_gev2_ms` Hessian plus explicit
   `R_wilson_from_background` rotation, scalar labels, and basis IDs.
   When full labels are supplied, all **four** quartic axes are permuted to
   Higgs-first canonical order with the mass-basis heavy-field ordering checked.

**Mandatory shared producer metadata** on both JSONs: `model`,
`matching_context_id`, `matched_seed_id`, `lagrangian_source_id`,
`mu_f_gev`, `mu_s_gev`, `scheme='MSbar'`, `units='GeV'`, and
`synthetic_fixture` (explicit boolean). The same context/seed/source ID must
hold across A and B. Wilson output must already include the evolved EFT1 C5,
and must explicitly exclude a separately added LLSS→C5 log. The full quartic
must declare `quartic_potential_convention='V4_lambda_real_symmetric_over_4factorial_v1'`
and `quartic_inventory_status='full_tensor'`; B must identify both scalar
basis IDs and supply the **orthogonal real** bridge `s_A=R s_B`.

**Important blocker:** As of the inspected main, the native B trajectory
`background.npz` does **not** itself contain the complete quartic tensor,
full mass Hessian, or scalar basis rotation needed by C. Passing that native
NPZ alone should fail as `BLOCKED`, not produce a fake T3 neutrino mass. The
new planned Stream 1/2 producer contracts can output the full tensors under
the exact paths above, and the same C command then consumes them without
modifications. No real A(-4) result is claimed until those producer files
actually exist and pass C preflight.

## Install and run (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Join-Path $env:USERPROFILE 'Downloads\t3_rough_stream3_threshold.zip'
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

# Check against the real local C/SM package and all synthetic fixtures.
.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_rough_stream3.py `
    tests\test_t3_parallel_c_threshold.py `
    tests\test_t3_parallel_c_real_integration.py

# First genuine A(-4) source-based case (never a fixture unless test opt-in).
.\.venv\Scripts\python.exe scripts\RunT3RoughAll16.py --models A_-4

# Then process the 16 models in one command, A(-4) attempted first.
.\.venv\Scripts\python.exe scripts\RunT3RoughAll16.py --overwrite

# View successes and actual blockers.
Import-Csv output\parallel_t3\rough_results\rough_fermion_first_v1\summary.csv |
    Format-Table model,status,stage,blocker,m1_eV,m2_eV,m3_eV -AutoSize
```

If the new Stream 1/2 scripts are not already run, this CLI can invoke their
*known actual* Python entrypoints, **only when input files are missing**, via
`--stream1-command '... --model {model} ...'` and
`--stream2-command '... --model {model} ...'`. No arbitrary default upstream
commands are assumed; no Matchete, RGBeta, or archived symbolic tasks run by
default. The command templates are executed without a shell in repository root.

For an explicitly conditional projection of an EWSB-resummed mass spectrum
onto the dimension-five Weinberg template, add `--allow-projection-residual`.
This does not make the SU(2)-noninvariant remainder physical SMEFT matching,
and the choice is recorded in the manifest. Default is to fail preflight.

A single C-only fixture smoke test (outputs **not** genuine):

```powershell
# Existing deterministic fixture generator, used only by pytest.
.\.venv\Scripts\python.exe -m pytest -q `
    tests\test_t3_rough_stream3.py -k 'synthetic_smoke_A4'
```

## Outputs

Under `output/parallel_t3/rough_results/rough_fermion_first_v1/`:

- `<slug>/C5_smeft_at_MS.npz`: EFT1 C5, finite scalar bubble, independently
  supplied hard C5 or **explicit zero approximation**, their sum.
- `<slug>/scalar_mass_basis.npz`: full Hessian, orthogonal rotation, all
  positive heavy scalar mass eigenvalues and eigenvectors.
- `<slug>/C5_smeft_low.npz`: evolved C5, complex Majorana matrix, sorted
  neutrino masses, Takagi and PMNS matrices when defined.
- `<slug>/observables.json`: sorted masses, Δm², mixing angles, J and δ if
  numerically defined, and precise basis conventions.
- `<slug>/manifest.json`: `RoughNumericalBenchmark` or explicitly synthetic
  smoke status, actual input SHA256s, EFT order, scales, C-only omissions,
  one-loop backend diagnostics, all finite bubble pairs and quartic channels.
- `summary.json` and `summary.csv`: all requested model results or actual
  blockages. Blocked models are recorded in `_blocked/<slug>.json` too.
- `_adapted_inputs/<slug>/`: **C-owned** normalised, checksum-stamped full
  JSON/NPZ input pair for reproducibility; never modifies upstream A/B.

Missing inputs produce no `C5_smeft*.npz` and no fabricated numbers. The runner
exits nonzero if any requested model is blocked; the CSV/JSON are still saved.
