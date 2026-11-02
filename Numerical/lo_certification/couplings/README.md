# T3 Workstream 2 — physical EFT1 coupling and scalar-basis certification

**Source inspected:** GitHub `Frankyang2002/T3neutrinomassmodelpipeline`, `main`, tree
`2c29e5ba2c0350b7e9be292e02ceb8bc7c52c0d3`, 2026-11-02, **read-only**.
The batch records the local Git blob hashes before producing source-backed outputs.
For the two observed drifting files only, it now performs a fail-closed *local*
structural re-audit instead of demanding byte-for-byte equality with GitHub:
`ScalarOnlyTensorAdapters.py` re-parses every selected native Matchete quartic
seed and compares **every symbolic channel and real-component derivative**
against the source-hashed v4 CG cache, while
`general_rgbeta_quartic_registry.json` validates all additional HHSS identities
against the original candidate inventory, forbidding removal of any older
HHSS candidate. New candidates are explicitly recorded as unnormalized.
All other drift remains a blocking error. These are checks for conditional
source compatibility, **not** a certificate for potential completeness or
Matchete/RGBeta normalization.
No external source or existing solver/output is modified. The implementation runs
against the user's real repository; the working checkout is not included in the ZIP.

## Fixed physics

`UV -> integrate out ALL heavy BSM fermions TOGETHER at mu_F -> EFT1 with
both heavy scalar multiplets (and all SM fermions) -> integrate out BOTH
heavy scalar multiplets TOGETHER at mu_S -> SMEFT`.

The tree-level Matchete `C_LLSS` is a **formal spinor polynomial**, not a
physically normalized Wilson coefficient. The signed `GammaCC -> Weyl`
factor and complete one-loop hard matching are NOT certified by the current
source. The source registry contains additional symmetry-allowed quartics
whose relative Matchete/RGBeta tensor normalizations have NOT been checked.
**Consequently no physical C5, physical neutrino mass or certified master
seven-term solver input is manufactured here.**

## Bundle files, repository-relative

- `Numerical/lo_certification/couplings/__init__.py`
- `Numerical/lo_certification/couplings/basis.py`
- `Numerical/lo_certification/couplings/batch.py`
- `Numerical/lo_certification/couplings/interfaces.py`
- `Numerical/lo_certification/couplings/README.md`
- `scripts/RunLOCouplingCertification.py`
- `tests/test_lo_couplings_basis.py`
- `tests/test_lo_couplings_batch.py`
- `tests/test_lo_couplings_integration.py`
- `output/general_observable/lo_certification_couplings/SOURCE_ONLY_PRECHECK_20261102.json`
- `output/general_observable/lo_certification_couplings/OUTPUT_SCHEMA_V1.json`

## PowerShell (from checkout root)

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WS2_EFT1_Coupling_Batch.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p 'test_lo_couplings_*.py' -v
python -m scripts.RunLOCouplingCertification --repo-root . --run-name ws2_all16_v1
Get-Content 'output\general_observable\lo_certification_couplings\ws2_all16_v1\summary.json' -Raw
```

**Do not** run `--allow-local-version-drift` as a physical certification: it
bypasses even the new local structural re-audit. With the normal/default run,
changed parser output is compared directly with hash-verified quartic caches,
and missing prior registry invariants cause a hard failure. Use a new
`--run-name` for each execution; output remains exclusive-write.
The batch will *not* overwrite existing output directories. Exit code 2 means
one or more real input source/interface stages were blocked; `summary.json`
and model-specific `failure.json` preserve that failure without substituting
zero couplings.

Optional pre-existing **fully enumerated diagnostic** inputs for the installed
seven-term + simultaneous scalar solver can be checked using the original
solver, without re-implementing its beta functions:

```powershell
python -m scripts.RunLOCouplingCertification --repo-root . `
  --run-name ws2_existing_solver_interface_v1 `
  --existing-master-coupling-root 'output\general_observable\lo_eft1_mixing\my_explicit_inputs' `
  --existing-joint-background-root 'output\general_observable\lo_joint_scalar_matching\my_backgrounds'
```

Those roots must be **existing input directories**, not the provisional
`conditional_numeric_tensors.npz` files from this workstream. The existing
`IntegratedEFT1JointScalarBridge.evaluate` demands `couplings.json` and
`<slug>.json` plus sibling NPZ with its original schemas, complete source
channels, SU(2)/U(1) gauge generators and Yukawa tensors. It returns a
conditional output **per unknown signed** `N_Gamma`, not a physical mass.

## Exact output semantics

For each model `A_m4` through `E_p2`:

| File | Description | Certification |
| --- | --- | --- |
| `symbolic_hhss.json` | Exact native HHSS components from existing Matchete-CG v4 cache, source SHA, all extra HHSS registry candidates | Source tensor verified, **incomplete general potential** if extra channels exist |
| `running_couplings_muF_muS.json` | All original renormalizable EFT1 running-state values at both scales, including SM Yukawas, gauge and quadratic parameters | Actual repository **illustrative** running, not certified physical parameter trajectory |
| `conditional_numeric_tensors.npz` | Native candidate quartics, `R_F`, `R_S`, relative rotation, full positive unbroken-EW illustrative Hessian at `mu_S`, source C tensor (formal) | **Conditional**, missing independent quartics not filled |
| `conditional_numeric_metadata.json` | All candidate names, masses, basis and scale assumptions, readiness blockers | Reject as physical input; deliberately **not** the existing joint-matcher schema |
| `existing_solver_contract.json` | Full existing integrated seven-term / joint-matching required input contract | Strict, source-derived missing-input contract |
| `status.json` | `verified/conditional/blocked` per stage, all physical blockers, file hashes and basis labels | Physical status **blocked** until proofs provided |

The global `source_versions.json` uses SHA-256 of local bytes and Git blob
SHA-1 to compare directly to the inspected main tree. `summary.json` has exactly
16 rows, including errors, with separate physical and provisional coupling
statuses. Only files inside `output/general_observable/lo_certification_couplings/`
are created during runs.

## HHSS normalization and scales

The exact SymPy source parser retained from
`RGE.running.intermediate.ScalarOnlyTensorAdapters.build_eft1_quartic_tensor`
produces the fourth derivative of the potential (not an arbitrary portal
coefficient):

\[
 V_4=\frac{1}{4!}\lambda_{abcd}\phi_a\phi_b\phi_c\phi_d,
 \quad \lambda_{abcd}=\frac{\partial^4 V_4}{\partial\phi_a\partial\phi_b\partial\phi_c\partial\phi_d}.
\]

Each independent Matchete coupling remains a named channel; numerical
coefficient tensors are **added** across channels without discarding suppressed
terms or source-level cancellations. All real Higgs legs and heavy scalar
legs are included (not just the neutral-pivot elements).

Let the source Matchete original real scalars be `s_0`. The fermion-boundary
coordinates at `mu_F` and instantaneous scalar basis at `mu_S` are

\[
 s_0=R_Fs_F=R_Ss_S,\qquad U_{F\leftarrow S}=R_F^T R_S,\qquad
 M^2_{F\,\mathrm{basis}}(\mu_S)=R_F^T M^2_0(\mu_S)R_F.
\]

The quartic and Wilson tensors are rotated on **both** heavy scalar slots
with the same orthogonal `R_F`. When `mu_F` and `mu_S` have different quadratic
mixing, `R_F` and `R_S` are constructed from their own endpoint inputs; the
mass matrix remains non-diagonal in the fixed fermion-boundary basis, and the
existing simultaneous bubble diagonalizes it once. The native six
self-bilinear `Minus/Plus` conventions and B/C mixed epsilon map are reused
from `Numerical/joint_scalar_matching/special_mass_v6.py` and are **explicitly
flagged unverified** when only the illustrative running state is available.
No `U_F = U_S` assumption is imposed.

For fixed overall one-loop order, coupling variations within an existing
one-loop LLSS-induced `C5` source contribute at one further loop. This does
not authorize replacing unverified `mu_F` values with `mu_S` ones in a claimed
physical renormalizable input. This batch records **both** states and leaves
the seven-term solver's full input format intact.

## What remains physically blocking

1. Independent signed Matchete `GammaCC` to master two-component Weyl
   normalization; complete hard all-fermion one-loop matching.
2. Real-component Matchete/RGBeta tensor equality (coefficients, signs,
   complex conjugate phases) for every candidate `lambdaH1Inv1/2`,
   `lambdaH2Inv1/2`, `lambdaT3` and all extra independent HHSS invariants.
3. Full RG-closed quartic and quadratic parameter sources in the existing
   RGBeta trajectory, including the extras flagged by the structural registry.
4. Complete SM-Weyl fermion basis and all SM Yukawa/gauge generator matrices
   in the **existing** `GeneralEFT1C5Mixing` normalizations; do not use its
   `synthetic` toy-gauge-zero mode as a physical prediction.
5. Independent basis phase / epsilon certification for mixed B/C and six
   self-bilinear models using actual Matchete native fields at `mu_F`.
6. General one-loop EFT1 operator-class closure and non-LLSS hard matching,
   if needed for a full leading-nonzero neutrino-mass result.

The numerical preview establishes only **which source-derived transformations
can be computed** from the current running snapshots. It must not be
interpreted as a complete one-loop physical EFT1 coupling handoff.
