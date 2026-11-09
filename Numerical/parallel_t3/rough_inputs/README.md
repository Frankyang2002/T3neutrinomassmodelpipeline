# T3 Stream 1 — provisional saved-B-to-A/C numerical input handoff

Run ID: `rough_fermion_first_v1`. Adds files only; does not replace existing B,
A or C modules and does not alter GitHub.

## Verified upstream contracts

Inspected repository `main` on 2026-11-10: `background/api.py`,
`background/contracts.py`, `background/fermion.py`, `background/handoff.py`,
`wilson/interchange.py`, `wilson/basis.py`, `wilson/coordinates.py`,
`wilson/transforms.py`, `threshold/README_C.md`, `threshold/models.py`,
`threshold/pipeline.py`, and B/C mass bridge v21–v30. B's committed
`B_all16_qualitative_v1/summary.json` records 16 successful candidate
renormalisable backgrounds and 16 formal F-tree boundaries. That is **not**
16 certified A/C endpoints. This batch consumes the local saved files, using
B's hash/continuity checks; it never redoes RGBeta, Matchete or the UV/EFT1 ODE.

The EFT order is **UV → integrate out ALL fermions at `MF=100000` GeV → EFT1
with S1 and S2 → integrate out BOTH scalars together at `MS=1000` GeV → SMEFT**.
This stream stops before the Wilson ODE and scalar threshold matching. Mass
eigenvalues can differ from `MS`: C will retain the mass-dependent logarithms.

## Windows PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Join-Path $env:USERPROFILE 'Downloads\t3_rough_stream1_inputs.zip'
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_rough_stream1_inputs.py

# A(-4) first. These are explicitly selected rough, not certified, inputs.
.\.venv\Scripts\python.exe -X utf8 scripts\RunT3RoughStream1Inputs.py `
  --model A_-4 --assume-formal-identity --assume-sm-weyl-yukawa `
  --gamma-factor 1 --hard-c5-zero-benchmark --reduced-potential

# A(-4) is deliberately first again; --overwrite permits its regeneration.
.\.venv\Scripts\python.exe -X utf8 scripts\RunT3RoughStream1Inputs.py `
  --all --overwrite --assume-formal-identity --assume-sm-weyl-yukawa `
  --gamma-factor 1 --hard-c5-zero-benchmark --reduced-potential

$s = Get-Content output\parallel_t3\rough_inputs\rough_fermion_first_v1\summary.json -Raw | ConvertFrom-Json
$s.models.PSObject.Properties | ForEach-Object {
  [pscustomobject]@{ Model=$_.Name; Status=$_.Value.status; Blocker=$_.Value.blocker }
} | Format-Table -AutoSize
```

For a nonzero hard Weinberg boundary, use `--hard-c5-root <relative_directory>`
with `<relative_directory>/<slug>/hard_C5.npz` containing `C5_hard[3,3]`, and
`hard_C5_manifest.json` containing the same-run UV state SHA, matching scale,
model and NPZ SHA256 as in Workstream B's README. Do **not** also use
`--hard-c5-zero-benchmark`. Absence of hard C5 does not prove it vanishes.

Omit `--reduced-potential` to *fail* on any nonzero unmapped native quartic.
Omit either identity/Yukawa approximation flag to fail closed. All original
saved values and their time trajectories are retained, even when their
*component contribution* is zeroed in the selected reduced scalar potential.

## Output contract

`output/parallel_t3/rough_inputs/rough_fermion_first_v1/`

- `a/<slug>/saved_native_trajectory.npz`: full saved UV/EFT1 source coordinate trajectories, including UV heavy fermions.
- `a/<slug>/background.npz`, `background.json`: **A `t3-parallel-a-input-v1`** with full `[N,45,45,ns]` Yukawas and `[N,ns,ns,ns,ns]` quartics, scales, gauge running, checksum, source hashes, basis labels. `ns=12` except E, for which `ns=16`.
- `a/<slug>/boundary.npz`, `boundary.json`: A native formal LLSS array and explicit source map/factor, `C5_hard`, **and** pre-expanded provisional canonical `C_full_provisional[45,45,ns,ns]` including the whole Weinberg polynomial. The existing A loader reads the formal input to reconstruct the tensor; it does **not** trust the pre-expanded tensor as independently certified.
- `a/<slug>/assumptions.json`: per-model signed factor, source identifiers, the native-nonzero values explicitly omitted by reduced potential at MF and MS, maximum absolute values over the whole trajectory, and explicit zero choices.
- `c/<slug>/background.json`, `background.npz`: **C `t3-parallel-c-background-input-v1`**, including `mass_squared_gev2_ms`, `R_wilson_from_background`, `quartic_real_full_ms`, `ye`, `yu`, `yd`, SM couplings and common context IDs.
- `summary.json`: actual model outcomes and blockers of this local run.

`quartic_inventory_status='full_tensor'` describes a **fully expanded rank-4
real tensor for the explicitly selected reduced potential**. The separate flags
`original_RGBeta_potential_fully_mapped=false` and
`reduced_potential_explicit_zeros=true` prevent confusion with the original
full RGBeta invariant inventory. This choice makes C's existing v1 loader
usable for a numerical reduced-potential benchmark. Do not describe the
result as the full original T3 potential.

## Construction / numerical assumptions

1. `C_LLSS_formal[3,3,2,2,n,n]` is reused directly from the saved same-run B
   Matchete flavour/source tensors. `GammaCC->Weyl = +1` is an explicitly
   selected **unverified** signed normalization. `heavy_source_from_canonical`
   is explicitly set to identity; notably B(-1)/C(-1) formal X/Y are **not**
   independently aligned with RGBeta S1/S2. C's `R_wilson_from_background=I`
   uses this same provisional frame, not a certified physical rotation.
2. The 45 SM left Weyl fields are ordered by the existing `ModelBasis`.
   The retained SM running `yu/yd/ye` matrices enter the conditional
   `Q u^c epsilon H`, `Q d^c H*`, `L e^c H*` polynomial under
   `L_Y=-Y_ija psi_i psi_j phi_a/2 +h.c.`; Yukawa signs/frame are not
   independently certified against all original RGBeta inputs.
3. Six basic quartics reproduce `T3RGBetaModel.wl` lines 383-457 exactly as
   a declared-potential convention: `lambdaH/2(HdagH)^2` etc. Adjoint
   invariants use the source's SU(2) generator contractions; `lambdaT3` uses
   the declared symmetric SU(2) doublet-to-triplet CG tensor in a
   **conditional** orthonormal triplet frame. Extra quartics not mapped by
   these explicit constructions, including `lambdaT3GeneralQ...` and selected
   B(-1)/C(-1) phase-sensitive terms, are set to zero *only when
   `--reduced-potential` is supplied*. Their actual **nonzero** B running
   values are not erased, and every dropped key/value is reported.
4. `threshold.models.from_parameters` constructs the heavy real scalar
   Hessian from B's actual MS quadratic inputs, including mixed B(-1)/C(-1)
   doublets and six self-bilinear splittings; nonzero complex phases are
   retained. The Hessian is unbroken-phase, not EWSB-resummed. Both scalars
   remain together to the selected `MS`.
5. Native B RGEs were integrated with the complete original scalar
   parameter inventory. The **projected reduced** rank-4 tensor is evaluated
   on that trajectory without reintegrating a changed reduced-potential RGE;
   therefore it is an intentionally approximate benchmark and not a
   self-consistent reduced model RG trajectory. No synthetic numerical
   parameters are substituted for saved B entries.
6. The hard Weinberg coefficient is **zero solely by explicit benchmark
   selection**, or an externally supplied symmetric nonzero `3x3` matrix.
   Neither option asserts the missing one-loop hard terms vanish. The scalar
   matching, any Wilson anomalous dimensions, and SMEFT evolution are out
   of scope.

`tests/test_t3_rough_stream1_inputs.py` exercises quartic normalization and
permutation symmetry across 16 model geometries, every-family SM Yukawa
embedding, mixing/splitting Hessians, omitted-term fail-closed behavior,
and the JSON/NPZ adapter I/O with *test-only* arrays. No real saved-B output
was executed in the sandbox without the original local `.npz` files.
