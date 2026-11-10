# A(-4) full-general source EFT1 reconstruction — v5

## Exact source evidence (GitHub `main` at `bad407a7`, 2026-11-10)

`output/parallel_t3/full_general_final_one_batch/thesis_a4_hermiticity_diagnosis_v4/A_-4/A_joint_scalar_lagrangian_audit.json` reports:

- Mass promotion reversible: `true` (MS1, MS2 from EFT order 1 to 0).
- `CheckLagrangian_original` and `CheckLagrangian_relabelled`: `False`.
- After removing W and lambdaT3 **for diagnostics only**, `CheckLagrangian` was still `False`.
- Old first-appearance alpha-key test claimed saved pure-SM and rebuilt UV LSM were not equivalent, but such a key is sensitive to Mathematica's dummy-dependent factor order.
- Log contains `CheckLagrangian::Hermiticity`. A physical match may NOT bypass this.

## What the full replacement does

New diagnostic path in `Lagrangian/RunT3AMinus4FullGeneralPhysicalScalar.wl`:

1. Rebuild SM, heavy scalar free Lagrangians and **all** scalar-only UV interaction monomials from `BuildT3Lagrangian`, after the existing heavy-field registry is restored.
2. Exclude UV terms containing NewFermion (already integrated out at mu_F).
3. Add **the two exact** matched F-tree LLSS operators with the independent flavour tensor W[p,q]. F hard O(hbar) remains separately saved and is not spuriously loop-matched again.
4. Compare complete rebuilt EFT1 with the exact matched-and-mass-promoted EFT1 using a canonical dummy-index graph key. For each SU(2), Lorentz or flavour index representation, the key minimises over permutations; it therefore is independent of dummy names and factor sort order while preserving all operator/coupling/CG content. Limit large searches instead of silently accepting equivalence.
5. Run CheckLagrangian on the reconstructed scalar base for diagnostics. Run CheckLagrangian on the **full reconstructed** EFT1 only if complete graph-equivalence is proven. Pass the rebuilt EFT to Match **only if both gates pass**.
6. Save diagnostics even if any gate fails, including term-multiset mismatches, so a subsequent repair can target a specific operator instead of toggling physics interactions.

No arbitrary scalar couplings or SM terms are dropped. This does not claim a physically finite C5: inherited F hard and EFT1 Wilson RGE are still absent.

## Files in ZIP (relative to project root)

- `Lagrangian/RunT3AMinus4FullGeneralPhysicalScalar.wl`: full replacement.
- `tests/test_t3_a4_source_reconstruction.py`: static checks.
- `docs/T3_A4_SOURCE_EFT1_CANONICAL_RECONSTRUCTION_v5.md`: these instructions.

## PowerShell (from project root)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_A4_SourceEFT1_Reconstruction_v5.zip" -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p test_t3_a4_source_reconstruction.py -v
python -X utf8 -m scripts.RunT3ThesisFinalOneBatchSafe --model A_-4 --run-id thesis_a4_source_reconstruction_v5 --debug
$dir = 'output\parallel_t3\full_general_final_one_batch\thesis_a4_source_reconstruction_v5\A_-4'
Get-Content "$dir\A_source_reconstruction_audit.json"
Get-Content "$dir\native_two_threshold_status.json"
```

## Interpretation

If `whole_EFT_graph_alpha_equivalent` and `complete_rebuilt_CheckLagrangian` are both true, the source-rebuilt EFT is used for one joint scalar match (S1 and S2). The overall run will **still** be `Blocked_CompletionBlocked` because this is scalar soft only. If the graph comparison fails, consult `mismatch.previews` to see the specific monomials retained, absent or differently expressed; *do not* override the equivalence gate.

If `scalar_base_CheckLagrangian` is already `False`, the underlying SM+scalar source or registry needs repair prior to any Wilson-specific matching. If base is `True` and full rebuilt is `False`, concentrate on the W/CG Hermiticity conventions. Until a complete UV→EFT1→SMEFT fixed-order matching and RG assembly is verified, treat all C5 outputs as provisional.
