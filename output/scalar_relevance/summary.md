# T3 scalar relevance — Batch 1 source audit

- Models: 16
- Extra quartic slots: 50
- Extra quadratic slots: 8
- Extra quartic slots observed in saved C5: 4
- Affected models: B_-1, C_-1

## Classification counts

| Classification | Rows |
|---|---:|
| higgs_vertex_not_seen_in_saved_C5 | 20 |
| leading_one_loop_topology | 16 |
| no_direct_single_vertex_LLSS_to_LLHH | 26 |
| not_direct_single_vertex_LLSS_to_LLHH | 48 |
| observed_in_saved_one_loop_C5 | 4 |
| one_loop_higgs_vertex_candidate | 32 |
| tree_level_scalar_spectrum | 32 |
| unsuppressed_tree_scalar_mass_effect | 8 |

## Caveats

- This is a graph-leg and saved-artifact classification, not a newly calculated C5 or a proof of generic zero couplings.
- Omitting a quartic from the specific direct one-loop LLSS->LLHH graph does NOT omit it from the full UV scalar potential or scalar RGEs.
- Higgs-bearing interactions not seen in saved matching are NOT proven irrelevant in another scheme, scalar basis or mass spectrum.
- Allowed quadratic coefficients affect tree-level scalar masses/mixing and cannot be dropped by one-loop power counting.
- One-loop RG generation followed by one-loop mass matching is formally higher order for zero-boundary couplings, but logarithmic enhancement and boundary-scale conventions must be checked.

## Input hashes

- `docs/general_scalar_trial_operator_manifest.json` SHA256 `25e07dec41a681d38b461272d136d41e2becd378b5bac4ded86832ea02d70582`
- `output/general_matching_symbolic_report.json` SHA256 `7042bf00134e9989943cbe978817103f8c56684ece7b879aa2d42ba988e63ec9`
- `output/general_scalar_audit_report.json` SHA256 `9a39288f1df7213857c07c993a1676b36f661299c8a2e76f527a8c85bbc3f04f`
- `output/general_t3_physics_summary.json` SHA256 `8cc68b6bf4a7f4e5ca2bae8083bbb512e293360026f701ce185e0e5471fa4f52`
- `output/general_t3_quadratic_eft1/quadratic_basis_audit.json` SHA256 `9a5b5a66ecec665c9b47093feddf2af2485738b979a6b21f6142f884df51d2e4`

See `operator_classification.csv` for a per-operator reason and omission criterion.
