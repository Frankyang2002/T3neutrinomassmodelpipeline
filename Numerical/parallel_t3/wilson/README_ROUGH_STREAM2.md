# Stream 2 — rough, full coupled EFT1 Wilson transport

**Only added files; no existing Wilson source is replaced.** Source baseline inspected via
GitHub connector, `Frankyang2002/T3neutrinomassmodelpipeline` public `main` at
`25eb5e13d09dd2f138176d4e36c4fc65edfa4bb6` (2026-11-09 UTC). Local uncommitted
files are not visible to this delivery. This uses the existing modules
`basis.py`, `coordinates.py`, `adm.py`, `interchange.py`, `workflow.py` and the
existing repository `RGE.running.intermediate.GeneralEFT1SixTermReferenceCandidate`.

## Physics and numerical ownership

- All fermions are integrated out **together** at MF=1e5 GeV in Stream 1.
- EFT1 evolves **both scalars** and a full `(45,45,ns,ns)` complex symmetric
  `C_psi2phi2` with `ns=4+2*dS1+2*dS2` in the same DOP853 solve.
- The beta is `16*pi^2*dC/dlnmu=beta(C,C*,Y,lambda,g)` evaluated from
  the existing candidate ADM and linearly interpolated full **running**
  canonical background in `ln(mu)`.
- The selected default `fonseca_eq22_six_term` is a **numerical candidate**,
  not independently certified. Alternative `ward_covariant_candidate` is
  opt-in. Known noncovariant `source_v1_unstarred` is disabled here.
- `C5` is never separately evolved in EFT1. Its full Weinberg polynomial is
  embedded in the same master tensor and projected via the existing
  `project_c5`. The convention is
  `C[nu_p,nu_q,H[2,R],H[2,R]] = C5[p,q]/4`.
- Stream 3 is responsible for integrating **all scalars together at MS**,
  matching to SMEFT, and evolving SMEFT. This export does **not** add any
  direct LLSS->C5 leading logarithm a second time.

## Canonical Stream 1 input, per model

`output/parallel_t3/rough_inputs/rough_fermion_first_v1/a/<slug>/`
where `<slug>` is `A_m4`, `B_m1`, `D_p2` etc. The adapter accepts either:

1. `canonical_input.npz` + `canonical_input.json` (also named
   `rough_input.npz` + `rough_input.json`) with all canonical background
   and Wilson keys in a single NPZ; or
2. `background.npz` + `background.json` (or `background_manifest.json`),
   and `boundary.npz` + `boundary.json` (or `boundary_manifest.json`).

The **full background**, never only invariant coupling names:

| Array | Shape, convention |
|---|---|
| `mu_gev` | `[N]`, endpoints exactly `1e3` and `1e5` GeV; either order |
| `Y` | `[N,45,45,ns]` complex, Weyl-symmetric |
| `lambda4` | `[N,ns,ns,ns,ns]` real and fully permutation-symmetric; `V4=lambda4 phi^4/4!` |
| `g2`, `gY`, `g3` | `[N]` real, in the documented non-GUT hypercharge convention |

The **tree-level fermion-first Wilson boundary** is exactly one of:

- `C_psi2phi2_tree_mf` or `C_full`: full complex `(45,45,ns,ns)`
  tree tensor; metadata `C_psi2phi2_tree_mf_is_canonical:true` or
  `C_full_is_canonical:true` is required.
- `C_LLSS_formal_muF`: complex `(3,3,2,2,ns-4,ns-4)` from a *declared*
  transformed formal Matchete boundary. It requires explicit
  `heavy_source_from_canonical` real orthogonal array `(ns-4,ns-4)`
  and `gamma_cc_to_weyl_signed_factor` in the boundary JSON.
  `--assume-formal-identity` explicitly selects the **unverified** identity/+1
  assumption for rough benchmarking only; it is recorded and not physical.

Optional hard Weinberg seed: `C5_EFT1_mf` **or** `C5_hard`, complex symmetric
`(3,3)`, is embedded once as the *whole* Weinberg polynomial. Do not provide
one of those keys if the full tensor already contains C5. If no C5 is
provided, `--c5-boundary zero` explicitly selects zero at MF. Without this
flag the adapter blocks missing C5 rather than inventing a zero boundary.

**Manifests** must carry the actual existing matching context, never IDs
synthesized by this workstream:

```json
{
  "model": "A_-4",
  "basis_id": "T3_RGBeta_reference_Weyl45_real_interleaved_v1",
  "fermion_labels": ["exact 45 fields returned by build_basis('A_-4')"],
  "scalar_labels": ["exact fields returned by build_basis('A_-4')"],
  "full_SM_three_generation_Weyl_basis": true,
  "matching_context_id": "the actual ID exported by Stream 1",
  "matched_seed_id": "the actual F seed ID exported by Stream 1",
  "lagrangian_source_id": "the actual immutable source ID from Stream 1",
  "scalar_basis_id": "the exact scalar basis ID used by Stream 1 and Stream 3 B",
  "mu_f_gev": 100000.0,
  "mu_s_gev": 1000.0,
  "scheme": "MSbar",
  "operator_basis": "psi2phi2_no_factorial",
  "quartic_potential_convention": "V4_lambda_real_symmetric_over_4factorial_v1"
}
```

The four exact IDs, labels, and scales may be in either of the two JSON
manifests, but mismatches block. SHA256 checks of the NPZ and optional
`source_sha256` of project-relative files are supported. An input with
`eft1__lambdaH`, `eft1__ye`, etc. **without** canonical `Y/lambda4`
blocks rather than incorrectly treating native invariant slots as a real
component tensor.

### Existing saved Workstream B (diagnostic path)

The old B data under `output/parallel_t3/background/B_all16_qualitative_v1/`
use the **different** contracts `t3-parallel-B-background-v1` and
`t3-parallel-B-fermion-boundary-v1`:
`background.npz` has `eft1__<RGBeta beta symbol>` trajectories, while
`fermion_boundary_formal.npz` has `C_LLSS_formal` in source Matchete GammaCC.
They do **not** contain complete canonical `Y/lambda4` or an independently
established scalar/spinor coordinate mapping. The `--probe-b-only` mode reads
these real saved records and checks source-state consistency, SHA256 and any
available nonzero formal LLSS norm. It cannot honestly produce a physical
C5(MS) from these arrays alone. Use Stream 1's matching canonical adapter
output to run the full Wilson ODE, without changing the native B archive.

## Outputs (only on a successful genuine canonical input run)

`output/parallel_t3/rough_wilson/rough_fermion_first_v1/<slug>/`:

- `wilson.json`: **exact** Workstream C `t3-parallel-c-wilson-input-v1`
  fields, inherited IDs, source hashes, ADM/backend, solver controls, norm
  trajectories, 7 named beta channel projections/norms and uncertainty flags.
- `wilson_at_MS.npz`: `C_psi2phi2_ms`, `C_psi2phi2_tree_mf`,
  `C5_EFT1_ms`, `C5_EFT1_mf`, `C_LLSS_ms`, `C_LLSS_mf`, full 16pi2 beta
  at both endpoints, per-channel C5 beta tensors and sampled norm arrays.
- `batch_summary.json` (top-level): disposition of the attempted models.

In `rge_improved` mode Workstream C consumes the full evolved tensor and
separately does scalar matching. The optional C `strict_one_loop` mode remains
**unavailable** because this package deliberately does not invent its
`C_psi2phi2_ms_strict_1loop` endpoint.

The final tensor is **not** a synthetic fixture. If any required genuine
canonical field is missing, there is a `blocked` record in the summary and
**no fake Wilson endpoint**. Source assumptions and known ADM limitations
are serialized to every successful model's manifest.

## Installation and execution (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

$zip = Join-Path $env:USERPROFILE 'Downloads\t3_rough_stream2_wilson.zip'
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

$python = if (Test-Path '.\.venv\Scripts\python.exe') { '.\.venv\Scripts\python.exe' } else { 'python' }
& $python -m pytest -q tests\test_t3_rough_stream2_wilson.py

# Diagnose existing genuine but native-only saved B records, without fabrication.
& $python -m scripts.RunT3RoughStream2 --all --probe-b-only

# Run A(-4) on the NEW canonical Stream1 output first, with explicitly zero C5(MF)
# only if that is the chosen boundary and no hard C5 was supplied:
& $python -m scripts.RunT3RoughStream2 --model A_-4 `
    --c5-boundary zero --adm fonseca_eq22_six_term `
    --rtol 1e-6 --atol 1e-21 --steps 6

# Attempt A(-4) followed by the remaining 15, preserving each real outcome.
& $python -m scripts.RunT3RoughStream2 --all `
    --c5-boundary zero --adm fonseca_eq22_six_term `
    --rtol 1e-6 --atol 1e-21 --steps 6 --overwrite
```

For an independently provided nonzero C5(MF) boundary, **omit**
`--c5-boundary zero` (the default is `require`). To rerun with changed
controls use `--overwrite`; this replaces only Stream 2-owned model outputs.

**Evidence from this delivery:** See the accompanying execution notes in
`STREAM2_EXECUTION.txt` in the archive, which distinguish isolated adapter
unit tests from the unexecuted Windows/native 16-model integration.
