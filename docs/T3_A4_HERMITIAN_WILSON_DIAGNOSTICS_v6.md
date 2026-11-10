# A(-4) EFT1 Hermiticity and source-preserving scalar probe (v6)

## Evidence and scope

The v5 local audit reported:

- `scalar_base_CheckLagrangian = True` after reconstructing the renormalisable scalar theory from the validated UV builder.
- `whole_EFT_graph_alpha_equivalent = false`; 10 mismatched monomial forms, including SU(2) CG orientation and scalar quartic identities.
- A mismatch in the SM Yukawa up-quark sector: `CG[eps[SU2L]]` vs `CG[Bar[eps[SU2L]]]`.
- The first threshold is the same F-only three-flavour integration. The full one-loop F hard source is preserved separately and must *not* be fed to a new one-loop scalar match.

**A mismatch of a CG head is not proof of physical inequivalence**; conversely, it is not a licence to identify two tensor contractions without checking their representation and numerical tensors. The existing exact graph certificate remains mandatory for authoritative matching.

## New checks

1. Always execute `CheckLagrangian` on `UV scalar base + saved W pair` (diagnostic only), even when graph equivalence fails.
2. Extract exactly one holomorphic `Proj[1]` Wilson term and one conjugate `Proj[-1]` term from the saved F boundary.
3. Apply Matchete's `PlusHc` to the **unchanged** holomorphic Wilson term. Require exactly one holomorphic and one antiholomorphic term, and exact alpha-graph equality of the holomorphic component with the saved F threshold.
4. Independently verify Hermiticity of `UV scalar base + PlusHc(matched holomorphic W)` using `CheckLagrangian`. Export this as a *diagnostic* WXF.
5. Record the syntactic CG orientations in the old SM, UV SM, old Wilson pair and generated pair.
6. Preserve the distinction between exact-boundary matching and a source-based exploratory match.

## Local commands

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Expand-Archive -LiteralPath "$HOME\Downloads\T3_A4_Hermitian_Wilson_Probe_v6.zip" -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p 'test_t3_a4_*.py' -v

# Audit-only, no provisional replacement when equivalence remains unresolved:
python -X utf8 -m scripts.RunT3ThesisFinalOneBatchSafe --model A_-4 --run-id thesis_a4_hermitian_wilson_audit_v6 --debug

# Optional one-off exploratory scalar-soft match in the source-built Hermitian theory:
$env:T3_A4_EXPLORATORY_SOURCE_MATCH = '1'
python -X utf8 -m scripts.RunT3ThesisFinalOneBatchSafe --model A_-4 --run-id thesis_a4_hermitian_source_probe_v6 --debug
Remove-Item Env:\T3_A4_EXPLORATORY_SOURCE_MATCH -ErrorAction SilentlyContinue
```

**Warning**: Exploratory mode can only start matching if the source scalar base checks, the original F holomorphic W is exactly retained, Matchete derives both chiralities, and the entire source Hermitian candidate passes CheckLagrangian. It does not establish equivalence to the saved F EFT1; its numerical output must not be used as physical C5. This option is deliberately off by default.

Inspect:

```powershell
$dir = 'output\parallel_t3\full_general_final_one_batch\thesis_a4_hermitian_wilson_audit_v6\A_-4'
Get-Content "$dir\A_source_reconstruction_audit.json"
Get-Content "$dir\native_two_threshold_status.json"
```

## Outcomes

- `reconstructed_full_CheckLagrangian_DIAGNOSTIC_ONLY = "True"`: the saved W pair is Hermitian when combined with the validated UV scalar base. The failed saved EFT check is then likely caused by representation/canonicalisation, not that Wilson pair alone.
- `... = "False"` but `source_Hermitian_W_complete_CheckLagrangian_DIAGNOSTIC_ONLY = "True"`: the saved W partner is suspect; the canonical Hermitian completion may provide exploratory scalar-soft information (not a certified boundary).
- Both fail: no scalar matching is attempted. The CG and chirality inventories distinguish the remaining issue.
- `W_PlusHc_pair_matches_saved_tree = true` and full graph equality still false: examine remaining SM and scalar CG mismatch. Full equivalence still requires a tensor-algebra proof, not only matching couplings.

## Unchanged scientific limitations

- All three F flavours are integrated together at `mu_F`; the independent heavy masses enter `W_12[p,q] = sum_r y1[p,r] y2[q,r]/MF[r]`.
- Both scalars are integrated together at `mu_S`.
- F hard one-loop operators are **not** yet inherited; EFT1 Wilson RGEs and finite renormalisation are **not** implemented; no physical full-general `C5` or Yukawa fit is certified.
- D(-2) is untouched.
