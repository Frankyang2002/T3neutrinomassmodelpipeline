# Workstream 4 — Native fermion/spinor evidence, batch 7

## Scope and epistemic status

Additive, read-only, fail-closed source verification. Three actual native
RGBeta model registrations and three native restricted Matchete action builds
are **required** before the all-16 candidate Weyl-tensor comparison runs.
The local Python unit tests use *synthetic evidence*; they do not mean the
live Matchete/RGBeta measurements have already succeeded.

**This batch does NOT certify the physical cross-package component dictionary.**
It measures native registered **internal** Yukawa/mass tensor coefficients and
field ordering. It does not determine the four-spinor GammaCC numeric phase,
the RGBeta epsilon/triplet component axes, or the relative signed source
normalisation. The real beta tensor factor `Sqrt[2]` is *not* a physical
Yukawa Feynman-rule coefficient. The RGBeta MF tensor factor is also recorded
but not assumed to be the physical mass coefficient.

Previously completed scalar portal/mixing and full one-loop checks remain
unchanged. No fitting, RGE, matching, Weinberg operator or SMEFT code is
modified. Six symmetry-allowed extra Yukawas remain zero in the restricted
benchmarks by assumption, not by an operator-completeness proof.

## Native witnesses

`RunLODirectUVNativeFermionRegistry.wl` launches in a **fresh** Wolfram kernel
with RGBeta and records **simultaneously** the two Yukawa registrations and
the fermion mass registration for:

- `A_p0`: SU(2) doublet Dirac `F,Fc`, mass invariant `epsilon`.
- `B_m1`: singlet Majorana `F,F`, mass invariant `singlet`.
- `C_m1`: triplet Majorana `F,F`, mass invariant `delta`.

Checks actual native registry's `Fields`, `Chirality`, `Invariant`, `Coupling`
and `CouplingBar`. Does not evaluate RGEs. `RunLODirectUVNativeFermionAction.wl`
launches a **separate** Wolfram kernel with Matchete, forces the four T3
experimental flags to `0` only within that child process, builds the
**restricted** original Lagrangians and exports `LFree` and `LInt` InputForm,
checking baseline Yukawa and MF presence. No `Match[]` is called.

## Conditional all-left-handed reconstruction

`native_fermion_action_bridge.py` reads original source, 16 frozen benchmarks,
`component_support_audit.json`, `matchete_generator_analysis.json`,
`cg_covariance_audit.json`, and the 16 frozen inputs. It does not import the
production builder's `construct` or the one-loop engine.

For the declared **candidate** basis, it independently reconstructs

- Dirac: symmetric off-diagonal `MF` flavour blocks in `[F,Fc]`.
- Majorana: `MF` flavour blocks times `U_F.T @ U_F`, where `U_F` rotates the
  native SU(2) generators into the chosen T3 weight basis.
- `y1_HC`: conjugate scalar source and `conjugate(y1 * CG1)` with the
  declared `F^c`/Majorana mapping.
- `y2_HC`: field scalar source and `conjugate(y2 * CG2)` with the declared
  `F`/Majorana charge-pairing map.

It checks the original benchmark SHA256 and both original CG/generator
provenance hashes for each model. The source audits are in
`output/general_observable/lo_direct_one_loop/source_audits/`, **not**
`Reports/output/full/multiplicity/`.

The majorana/Dirac model coverage is 14 Dirac plus 2 Majorana. The latter are
`B_-1` and `C_-1`.

Successful tensor equality is proof **only within these explicitly declared
candidate Weyl conventions**. A numerical rephasing or an internal RGBeta
prefactor is not an independent 4-spinor action normalization certificate.
Certification stays `0/16`.

## PowerShell commands (repository root)

Install the ZIP into the repository and run:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
.\.venv\Scripts\python.exe -m pytest tests\test_lo_direct_uv_native_fermion_action_bridge.py -q

$native = 'output\general_observable\lo_certification_direct_uv\native_fermion_action_01'
foreach ($model in @('A_p0','B_m1','C_m1')) {
    wolframscript -file scripts\RunLODirectUVNativeFermionRegistry.wl `
        $model "$native\${model}_rgbeta.json"
    if ($LASTEXITCODE -ne 0) { throw "RGBeta native witness blocked: $model" }
    wolframscript -file scripts\RunLODirectUVNativeFermionAction.wl `
        $model "$native\${model}_matchete.json"
    if ($LASTEXITCODE -ne 0) { throw "Matchete native witness blocked: $model" }
}

.\.venv\Scripts\python.exe scripts\RunLODirectUVFermionSpinorBridge.py `
    --project-root . `
    --native-dir $native `
    --output-dir output\general_observable\lo_certification_direct_uv\native_fermion_spinor_bridge_01
```

Output:
`output/general_observable/lo_certification_direct_uv/native_fermion_spinor_bridge_01/native_fermion_spinor_bridge.json`.
Each native file is immutable: do not rerun a successful witness to the same
path; choose a new directory if necessary. The Python runner also refuses to
overwrite a nonempty output directory.

## Remaining blockers

1. Native numerical `GammaCC` orientation and phase relative to the declared
   four-component chiral basis.
2. Native RGBeta numeric SU(2) epsilon/triplet basis in the exact same axes.
3. Signed, complex-conjugation-aware **physical** `y1`, `y2`, `MF` coupling
   dictionary between the packages, including the Majorana `1/2` convention.
4. Relative `lambdaT3` coupling normalisation across the packages.

Correct EFT order remains:
**UV -> integrate out ALL heavy fermions at muF -> EFT1 with SM fermions and
BOTH BSM scalars -> integrate out both scalars at muS -> SMEFT**.
This batch is a broken-phase UV verification, not threshold matching.
