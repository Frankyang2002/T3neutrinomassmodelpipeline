# Existing 16-model fermion-first boundary → joint scalar matching

## Scope

New, isolated files: `Numerical/joint_scalar_matching/fermion_boundary_bridge.py`,
`scripts/BridgeT3JointScalarEFT1.py`, `tests/test_t3_joint_scalar_fermion_boundary_bridge.py`,
and this note. They do **not** alter upstream fermion matching, EFT1 mixing,
first-threshold modules, historical results, or GitHub. Only the established
`output/general_observable/lo_joint_scalar_matching/` subtree is writable.

Verified read-only GitHub `main` source commit as of 2026-11-02:
`e5e4a32f6f5a96eb92ac4fc2c9884127e56a2ce4`.

Existing all-16 release:
`output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/`

- 16/16 tree-generated, **source-resolved formal LLSS tensors at muF=100000 GeV**.
- All independent native PL operators, including small terms, are retained.
- PR Hermitian conjugates are not added again.
- Canonical row index is `[lepton SU2 i, flavour p, lepton SU2 j, flavour q, scalar a, scalar b]`, **1-based**.
- Real fields 1–4 are Higgs, heavy scalars start at index 5. The bridge maps
  to `C_LLSS[p,q,i,j,a-5,b-5]`, expanding the distinct lepton/scalar exchange
  orbit **without applying its size as another multiplicity**.
- Both B_-1/C_-1 mixed-doublet `X/Y` bases and six split-real cases are read
  from the actual explicit scalar field offsets. No implicit conversion to an
  S1/S2 basis or invention of an SU(2) tensor is permitted.
- Independently check that the sum over *every* PL source reproduces the saved
  source-summed tensor. Export that exact formal tensor plus every source term.

**Critical limitation:** The released `GammaCC → master Weyl` absolute sign and
factor are **not certified**; the record explicitly sets `certified=false` and
factor/sign `null`. Its `muF` tree tensor has **not been run to `muS`**. The
bridge never labels this provisional data a physical Wilson coefficient.

## PowerShell: ready now (read-only upstream; formal source audit)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
python -m unittest discover -s tests -p 'test_t3_joint_scalar_fermion_boundary_bridge.py' -v
python -m unittest discover -s tests -p 'test_t3_joint_scalar_matching.py' -v
python -m scripts.BridgeT3JointScalarEFT1 `
  --release-dir 'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f' `
  --output-dir 'output/general_observable/lo_joint_scalar_matching/fermion_boundary_audit_v1'
```

Expected: source-resolved, **not physically normalized**, `muF` tensors in
`fermion_boundary_audit_v1/muF_formal_sources/{model}.npz` and metadata JSON;
summary has 16 rows and `physical_C5_certified=false`. It intentionally
produces **no C5** because no complete normalized muS endpoint is present in
this handoff. Use a fresh output directory for reruns: old summaries are not
overwritten.

## Independent EFT1 muS endpoint contract (future input)

To run the physical-input **candidate bubble**, provide ONE JSON and ONE NPZ per
model, filenames `A_m4.json`, `A_m4.npz`, etc, in a separate read-only folder.
For a selected model `--model B_-1` only its pair is required.

Metadata JSON (no values may be inferred):

```json
{
  "format": "t3-joint-scalar-eft1-endpoint-v1",
  "model": "B_-1",
  "source_boundary_sha256": "<sha256 of B_m1/boundary.json>",
  "source_scalar_real_labels": ["<each of the 8 source real field labels in index order>"],
  "tensor_npz": "B_m1.npz",
  "mu_f_gev": 100000.0,
  "mu_s_gev": 1000.0,
  "reference_scale_gev": 1000.0,
  "renormalization_scheme": "MSbar",
  "tensor_convention": "T3_REAL_LLSS_HHSS_WEINBERG_v1",
  "scalar_mass_origin": "unbroken_electroweak",
  "quartic_channel_names": ["native_portal_1", "native_portal_2"],
  "eft1_evolution": {
    "status": "renormalized_full_EFT1_to_muS",
    "same_muS_wilson_quartic_mass_C5": true,
    "no_separate_scalar_threshold": true,
    "C_LLSS_used_in_bubble_order": "hbar0",
    "C5_EFT1_includes_running_once": true
  },
  "normalization_certificate": {
    "status": "independent_signed_amplitude_verified",
    "source_boundary_sha256": "<sha256 of B_m1/boundary.json>",
    "master_operator_convention": "T3_REAL_LLSS_HHSS_WEINBERG_v1",
    "gamma_cc_to_master_weyl_signed_factor": "<nonzero REAL NUMBER, not a string>",
    "evidence_file": "B_m1_signed_amplitude_evidence.txt",
    "evidence_sha256": "<sha256 of this independently produced file>"
  }
}
```

Every model must have independent signed-amplitude evidence, actually in the
same directory and hash-matched. **A certificate is a provenance assertion by
its supplier, not a mathematical proof supplied by the adapter.** The currently
published Matchete `GammaCC` probe reported failure to load the definition and
**does not satisfy this gate**. Do not fill it with a guessed `+1` or `-1`.

Required NPZ arrays:

| Name | Shape | Meaning |
|---|---|---|
| `source_to_matching_real_U` | `(n,n)` real orthogonal | `s_source = U @ s_matching` (required, including B/C and split-real cases) |
| `C_LLSS_muF_matching` | `(3,3,2,2,n,n)` complex | Signed source boundary at muF, rotated using U, **checked by this adapter** |
| `C_LLSS_muS_hbar0` | `(3,3,2,2,n,n)` complex | Tree-order Wilson tensor transported to common muS; no one-loop self-running inserted into an already one-loop scalar bubble |
| `lambda_HHSS_muS` | `(4,4,n,n)` real | All native quartic portal terms at same muS, fully normalized |
| `lambda_HHSS_channels_muS` | `(k,4,4,n,n)` real | Named, independent contributions, sum **exactly** to total quartic tensor |
| `mass_squared_gev2_muS` | `(n,n)` real symmetric | **Both** scalar multiplets together, all masses/bilinears/mixing |
| `C5_EFT1_muS` | `(3,3)` complex symmetric | Inherited direct Weinberg operator, EFT1 RG logs already included **once** |

`n = 2*dS1 + 2*dS2`: 8 for A/B/C/D, 12 for E. For a mixed
source `B_-1/C_-1`, the released basis is already `X/Y`; the endpoint author
must supply an explicit orthogonal `U` into the **target** real S1/S2
matching basis, with CG/epsilon signs from the same source conventions. For
split-real cases, `U` includes any required component permutation/rotation.
Identity is valid **only if independently justified** by matching field names
and component conventions.

The external endpoint must certify its own one-loop renormalized EFT1
transport. Existing historical symbolic expressions with literal `MF`, `MS`
and one-generation `y1/y2` are not this numeric tensor endpoint. No RG or
matching scale is inferred from those historical files.

Numerically the bridge checks the source-to-target rotation, the independent
source-reconstructed `muF` tensor, all tensor shapes/symmetries, scale, scheme,
positive eigenmasses and quartic channel completeness, then writes the
adapter's standard JSON/NPZ and calls `RunT3JointScalarMatching.compute`.
The resulting `C5` remains marked **`physical_full_one_loop_certified=false`**:
independent complete finite threshold counterterms/wavefunction or other EFT1
operators have not been proved by this bridge.

When endpoint files genuinely exist:

```powershell
python -m scripts.BridgeT3JointScalarEFT1 `
  --release-dir 'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f' `
  --endpoint-dir 'output/general_observable/lo_joint_scalar_matching/physical_muS_endpoints' `
  --output-dir 'output/general_observable/lo_joint_scalar_matching/conditional_muS_v1'
```

`--endpoint-dir` is a **read-only** input. If a model is absent, its certificate
is open, or numerical parity fails, no matching input/C5 is created for it and
that model is `blocked` in `summary.json`. Command exits nonzero for any block.

## Physics accounted for

The scalar bubble is *one* common `muS` threshold with all masses/eigenstate
pairs. All `muF` → `muS` logs belong in the separately supplied EFT1 `C5_EFT1`
and `C_LLSS` transport, **not** repeated as a second contribution to matching.
`C_LLSS_muS_hbar0` implements strict `O(hbar)` bubble counting. The single
bubble's finite `B0bar` and its explicit matching-scale logarithm are separately
reported, but are components of the *same* correction.

**Unresolved independently:** absolute four-component-to-two-component
spinor normalization; complete model-specific EFT1 RG evolution; independent
full one-loop SMEFT matching (other diagrams, finite field redefinitions and
counterterms). This bridge cannot certify a physical all-16 neutrino mass
spectrum without that upstream evidence.
