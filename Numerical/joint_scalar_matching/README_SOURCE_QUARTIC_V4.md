# v4 — Matchete-native HHSS tensor extraction and scalar endpoint assembly

**Project paths (new files only):**

- `Numerical/joint_scalar_matching/source_quartic_v4.py`
- `Numerical/joint_scalar_matching/README_SOURCE_QUARTIC_V4.md`
- `scripts/RunT3ScalarQuarticSourceBatchV4.py`
- `tests/test_t3_joint_scalar_source_quartic_v4.py`

This is the **first automated native scalar-potential tensor source extractor** in this worker. It consumes existing `output/full16_symbolic/T3_<class>_alpha_<tag>/data/eft1_after_F_scalar_quartic_seed.json` for all 16 models, with the existing real-tensor/CG parser `RGE.running.intermediate.ScalarOnlyTensorAdapters.build_eft1_quartic_tensor`. It does not invent a quartic, a scalar mass, an RG trajectory, or a Matchete/Weyl normalization sign.

## Source-only batch — usable immediately

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
python -m unittest discover -s tests -p 'test_t3_joint_scalar_source_quartic_v4.py' -v
python -m scripts.RunT3ScalarQuarticSourceBatchV4 `
  --output-dir 'output\general_observable\lo_joint_scalar_matching\source_quartic_v4'
```

The runner produces **actual source-derived symbolic tensors**, not 16 fake numeric quartics: `native_quartic_channels/A_m4.json` through `E_p2.json` plus `summary.json`. A channel is an independent Matchete quartic coupling, with complete complex-conjugate contributions and exact SU(2) CG coefficients, all factorial-normalized to

\[
 V_4=\frac{1}{4!}\lambda_{abcd}\phi_a\phi_b\phi_c\phi_d.
\]

Only coefficients with **two Higgs real legs and two heavy-scalar real legs** are exported. All independent quartic channels are retained, including those not coupled by the usual `lambdaT3` vertex. The parser transforms Matchete `L4` to `V4=-L4` before extracting coefficients. Each symbolic tensor is in the *original* `NewScalar1/NewScalar2` real basis. These are **not** yet numerical coefficients at `mu_S`, and the source quartic seed by itself does not prove completeness of the full general potential.

`source_basis_matches_fermion_boundary` distinguishes the eight ordinary cases from the six split-real plus B/C mixed cases. For split/mixed models the source quartic basis is **not** the fermion-boundary basis; the exporter does not silently identify them.

## Optional actual numerical muS endpoint

Supply one JSON per model in a separate directory, e.g. `scalar_endpoints_muS/A_m4.json`. These must be **actual same-scale** scalar quantities from a verified trajectory; UV `base_state` values at `1e7 GeV` cannot be copied to `mu_S` and relabeled as a physical EFT1 endpoint.

```json
{
  "format": "t3-source-quartic-scalar-endpoint-v4",
  "model": "A_-4",
  "seed_quartic_sha256": "COPY_FROM_SOURCE_QUARTIC_V4_OUTPUT",
  "source_boundary_sha256": "COPY_FROM_SOURCE_QUARTIC_V4_OUTPUT",
  "renormalization_scheme": "MSbar",
  "mu_f_gev": 100000.0,
  "mu_s_gev": 1000.0,
  "scalar_mass_origin": "unbroken_electroweak",
  "both_scalar_multiplets_active_until_single_muS": true,
  "masses_and_quartics_at_same_muS": true,
  "hhss_source_operator_set_reviewed": true,
  "quartic_couplings_muS": {
    "lambdaH1Inv1": 0.0,
    "lambdaT3": {"real": 0.0, "imag": 0.0}
  },
  "mass_squared_gev2_original_real_basis_muS": "REPLACE_WITH_REAL_SYMMETRIC_NUMERIC_N_BY_N_MATRIX",
  "provenance": {
    "scalar_coupling_trajectory": "REQUIRED_ACTUAL_SOURCE",
    "scalar_mass_trajectory": "REQUIRED_ACTUAL_SOURCE",
    "hhss_source_review": "REQUIRED_OPERATOR_INVENTORY_REVIEW"
  }
}
```

The example is **illustrative only, deliberately not a valid input**. The coupling dictionary must list **every coupling** in `couplings_required_at_muS` for the chosen model; explicitly zero couplings are accepted, but omissions are errors. The original-basis mass matrix must have the exact dimension `(8,8)` for A–D or `(12,12)` for E. A positive real spectrum is required.

For B/C mixed or split-real fields, also supply the checked orthogonal matrix `original_real_equals_boundary_real_R` with

\[
 s_{\rm original}=R\,s_{\rm fermion\ boundary},\qquad R^T R=I,
\]

and `special_scalar_basis_transform_externally_validated=true` with `provenance.original_to_boundary_basis_map`. This must include epsilon conventions and complex-to-real phases from the *actual* source. The exporter transforms **both** `M_S^2` and every quartic channel to the fermion boundary coordinates, and passes identity as `source_to_matching_real_U` to avoid double rotation. Ordinary models with identical real component layouts use identity automatically.

Once valid endpoint files exist:

```powershell
python -m scripts.RunT3ScalarQuarticSourceBatchV4 `
  --endpoint-dir 'output\general_observable\lo_joint_scalar_matching\scalar_endpoints_muS' `
  --output-dir 'output\general_observable\lo_joint_scalar_matching\source_quartic_muS_v4' `
  --run-v3
```

This writes `numeric_backgrounds/<slug>.json/.npz` in the **existing v2/v3 conditional matching input format**, then optionally invokes the installed `scripts.RunT3JointScalarBatchV3` to calculate the **per-unit-unknown-`N_GammaCC`** one-loop LLSS response, *not* a certified physical `C5`. The optional run fails closed if any requested endpoint is absent or invalid. Use `--model A_-4` for incremental numerical tests and choose a fresh output directory for each run; existing summaries and numeric backgrounds are not overwritten.

## Important limitations

1. The saved original-field quartic seed is a source of exact CG/potential **tensor components**, not proof that all symmetry-allowed general quartic interactions were included. The endpoint requires explicit operator-inventory review.
2. Existing original-field quartic seeds are **not** the physical split-real or mixed X/Y mass-basis quartics. There is no automatic extrapolation of experimental mass-basis data into a certified source-to-boundary rotation.
3. The signed native Matchete Feynman diagnostic has **not** established the absolute `GammaCC`-to-Weyl normalization or the complete hard one-loop F matching. v4 never guesses this factor.
4. The source-only run resolves a concrete missing **analytic quartic-CG tensor generation** step. Numerical source-matched scalar masses and running quartic values at `mu_S` remain upstream numerical dependencies.
5. EFT hierarchy: **UV → all F together at `mu_F` → EFT1(S1,S2) → all scalars together at one `mu_S` → SMEFT**. No sequential scalar thresholds.
