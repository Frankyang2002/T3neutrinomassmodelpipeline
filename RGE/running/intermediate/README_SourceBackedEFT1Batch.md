# EFT1 source-backed large implementation batch — no physical normalization assumed

**Source inspected:** live repository `main` on 2 November 2026, commit `73ce0f025ac14b8b5514d4f68098e149d8f48d26`. It contains the 16-model formal fermion-first release, the original seven-term `GeneralEFT1C5Mixing.py` and full-SU2 kernel. The separate joint scalar matching workstream already implements conditional finite bubbles; this batch DOES NOT duplicate or modify that code.

## New files and locations

- `RGE/running/intermediate/SourceBackedEFT1Batch.py` — strict actual-source-to-master adapter; all seven one-loop anomalous-dimension terms, source and quartic channel decompositions, full Weinberg component audit; conditional only.
- `scripts/RunSourceBackedEFT1Batch.py` — all-16 CLI with preflight and exclusive writes.
- `tests/test_source_backed_eft1_batch.py` — real-component layouts, complex flavour, source sums, seven-term and SU2 algebra, fail-closed input validation, 16-model CLI/no-overwrite tests.
- `RGE/running/intermediate/README_SourceBackedEFT1Batch.md` — this guide.

Requires **existing** files (from previous delivered and installed ZIPs):

- `RGE/running/intermediate/GeneralEFT1C5Mixing.py`
- `RGE/running/intermediate/FullSU2SourceResponse.py`
- `RGE/running/intermediate/ProvisionalSourceKernel.py`

These are imported read-only. No previous files or outputs are replaced. Source release remains unchanged.

## Ready-now command: generate all 16 genuine-source missing-input contracts

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_EFT1_SourceBacked_LargeBatch.zip" -DestinationPath .
python -m unittest discover -s tests -p 'test_source_backed_eft1_batch.py' -v
python scripts\RunSourceBackedEFT1Batch.py `
  --release-root 'output\general_observable\lo_fermion_boundary\releases\source_021dcc1e8fe8ba0f' `
  --stage22-diagnostics 'output\general_t3_stage22\stage22_aligned_geometry_diagnostics.json' `
  --output-root 'output\general_observable\lo_eft1_mixing\source_backed_batch_v1'
Get-Content 'output\general_observable\lo_eft1_mixing\source_backed_batch_v1\summary.json' -Raw
```

Expected: 16 `required_couplings_contract.json` files and one `summary.json`, `physical_C5_computed=false`. The optional Stage22 reader audits all 16 model entries and their extra-quartic candidate geometry while refusing to turn any candidate into numerical couplings. The recorded 46 blocks (50 symbolic slots) are explicitly *not certified physical quartic coefficients*.

This run has **real published Chat-1 LLSS sources**, not synthetic Wilson coefficients. It does **not** invent missing quartic/Yukawa/gauge values or a GammaCC factor, and therefore no numerical physical anomalous dimensions are output yet.

Output is restricted under `output/general_observable/lo_eft1_mixing/<new-run-name>/` and uses exclusive file creation; for reruns use a **fresh** output path. `--check-only` performs all input audits without writing. `--model B_-1` processes one model.

## Optional complete diagnostic coupling inputs (all seven terms)

With external couplings already calculated for **each selected model**, put `couplings.json` at `<coupling-root>/<tag>/couplings.json`. The exact required shape is specified by each generated `required_couplings_contract.json`; the following is a *schema sketch*, **not certified physics values**:

```json
{
  "format": "t3-eft1-explicit-master-couplings-diagnostic-v1",
  "model": "A_-4",
  "source_boundary_sha256": "EXACT_SOURCE_SHA256_FROM_CONTRACT",
  "mu_F_gev": 100000.0,
  "mu_S_gev": 1000.0,
  "n_fermions": 6,
  "n_real_scalars": 12,
  "renormalization_scheme": "MSbar",
  "master_conventions": {
    "master_wilson": "C_ijab psi_i psi_j phi_a phi_b + hc",
    "weinberg": "(1/2) C5_ij (L_i.H)(L_j.H) + hc",
    "quartic": "V4=lambda_abcd phi_a phi_b phi_c phi_d/4!"
  },
  "coverage_assertion": {
    "authority": "diagnostic_only",
    "all_active_weyl_fermion_components_enumerated": true,
    "all_yukawa_components_enumerated": true,
    "all_quartic_channels_enumerated": true,
    "all_gauge_generators_enumerated_or_explicitly_toy_zero": true,
    "toy_gauge_couplings_set_to_zero": true
  },
  "yukawa_components": [],
  "quartic_channels": [
    {"label": "EXAMPLE_UNPHYSICAL", "components": [{"indices": [3,3,5,7], "value": 0.1}]}
  ],
  "gauge_sectors": []
}
```

The example deliberately sets **toy** gauge interactions to zero and has an incomplete quartic dictionary, so it is not valid as a physical T3 point. For nonzero gauge couplings `gauge_sectors` must list named `SU2L` (three generators) and `U1Y` (one generator), with complete complex `fermion` and `scalar` matrices per generator; the existing master RGE checks matrix shapes and Hermiticity. The spinor normalization and certified full operator inventory remain unresolved even for a complete diagnostic coupling file.

```powershell
python scripts\RunSourceBackedEFT1Batch.py `
  --release-root 'output\general_observable\lo_fermion_boundary\releases\source_021dcc1e8fe8ba0f' `
  --coupling-root 'output\general_observable\lo_eft1_mixing\external_diagnostic_couplings' `
  --output-root 'output\general_observable\lo_eft1_mixing\source_backed_seven_term_v1'
```

The optional result is per-unit-unknown `N_Gamma`, including every PL source, *all seven* master-RGE terms, labelled quartic-source inventory, complex 3x3 flavour matrices, 21-point fixed-order running logs, independent full-SU2 quartic projection, and scalar-matching handoff metadata. It is **not** a physical C5 coefficient. A full-Ward mismatch is reported without promoting the neutrino-component pivot to an invariant coefficient. `source_to_master_fidelity` verifies every source against independent coordinate expansion and no PR source is added. Two independent quartic source channels with the same component **add**, never overwrite.

## Order counting / outstanding physics

The only matching sequence remains `UV -> integrate all fermions at muF -> EFT1(H,S1,S2) -> integrate both scalars together at muS -> SMEFT`. This module implements conditional **one-loop fixed-order EFT1 C5 anomalous-dimension logs only**; it does not perform the simultaneous finite scalar threshold (already a separate workstream), the hard-fermion one-loop matching, dim-5 operator-class closure or post-threshold SMEFT running. The actual signed `Matchete GammaCC -> master Weyl` factor is **unknown** (`null`). The new `matchete_gamma_cc_evidence_v2.txt` is source evidence, not a signed amplitude calibration. Stage 22 quartic geometry explicitly says `physical_coupling_map_certified=false`, so it is not silently treated as a numerical quartic dictionary.

At one-loop order for a tree LLSS source and one-loop scalar bubble, LLSS one-loop self-running inside that bubble enters at two loops; do not count it again. Conditional fixed-coupling EFT1 logs scale as `ln(muS/muF)/(16*pi^2)`.
