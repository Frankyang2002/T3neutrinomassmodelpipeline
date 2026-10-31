# Workstream C — native Y-threshold SU(2) coefficient comparison (experimental)

## Purpose

This audit advances B(-1)/C(-1) from coupling-monomial support to a source-level **component-coefficient comparison**. Its inputs are *already computed* EFT1 LLSS, EFT1 HHSS and Matchete Y/LLHH hard-pole expressions. No production matching files are replaced.

This audit deliberately has **two stages** because the CG tensors currently used in the local Y-matching source are Matchete kernel registrations, not assumed 2x2 epsilon/delta arrays:

1. **Native CG export** uses a modified, standalone **copy** of the existing Y-only experimental source, with an early-exit hook before `RunT3Matching`. It reuses the original Lagrangian builder, original mass-basis input and the same independent canonical Y-only setup. No `Match` call is performed in this invocation. The original script is never edited.
2. **Exact Python contraction** interprets the exported Matchete CG arrays with exact SymPy numbers, evaluates external `l_SU2=1, H_SU2=2` components, sums the internal X/Y SU(2) indices, performs complex-scalar Wick pairings with conjugation, and compares monomial coefficients to the saved raw Y/LLHH pole residue (not to a renormalized Wilson coefficient).

The common `1/M_F` factor and Matchete `hbar` flag are stripped for coefficient comparison. This is **not** a derivation of an MS-bar counterterm; those factors and the Majorana bilinear normalization must be handled independently in final physical matching. A ratio between pole and tree contraction is a *diagnostic*, not a renormalization prescription.

## Source files consumed (local project paths)

- `Lagrangian/RunT3BCYOnlyCanonicalReconstruction.experimental.wl` (read only; inspected blob `57afcd888598f86cccc02c069916fee694f0aed3`)
- `output/general_t3_mass_eigen_matchete/{B,C}_m1_final_parameters.json`
- `output/general_t3_mass_eigen_matchete/eft1_fermion_first_tree/{B,C}_m1_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt`
- `output/general_t3_mass_eigen_matchete/eft1_fermion_first_tree/{B,C}_m1_final_eft1_fermion_first_tree_llss_tree.inputform.txt`
- `output/audits/t3_bc_y_renormalization_ledger_v1/{B,C}_m1_y_llhh_raw_pole_residue.inputform.txt`
- Previous v2 audit script: `scripts/AuditT3BCYMixedScalarPoleSupport.experimental.py`

## Files installed by this ZIP

- `scripts/ExportT3BCYNativeCG.experimental.py` — prepares a fresh-kernel **copy** of current local Wolfram source with a native CG export hook.
- `scripts/CompareT3BCYNativeSU2Coefficients.experimental.py` — exact B/C SU(2) contractions and pole-ratio diagnostic.
- `docs/WorkstreamC_YNativeSU2Component_README.md` — this guide.

## PowerShell (repository root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

python scripts\ExportT3BCYNativeCG.experimental.py --self-test
python scripts\CompareT3BCYNativeSU2Coefficients.experimental.py --self-test

python scripts\ExportT3BCYNativeCG.experimental.py --prepare

$env:T3_BC_Y_INDEPENDENT_CANONICAL = "1"
Remove-Item Env:T3_BC_Y_INDEPENDENT_CHECK_ONLY -ErrorAction SilentlyContinue

$dir = "output\audits\t3_bc_y_native_cg_components_v1"
$base = "output\general_t3_mass_eigen_matchete"

foreach ($stem in @("B_m1", "C_m1")) {
    wolframscript -file "$dir\${stem}_ExportNativeCG_NO_MATCH.experimental.wl" `
        "$base\${stem}_final_parameters.json" `
        "$base\eft1_fermion_first_tree\${stem}_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt" `
        "$dir\${stem}_native_CG.json"
    if ($LASTEXITCODE -ne 0) { throw "Native CG export failed: $stem" }
}

python scripts\CompareT3BCYNativeSU2Coefficients.experimental.py `
    --output "$dir\comparison.json"
```

The Python preparation deliberately fails if the original Wolfram source has a different SHA from the inspected Git blob; check the current source before using `--allow-changed-source`. This is not a production source modification.

## Status and limitations

The supplied script self-tests have been run in a Python/SymPy environment. **A Wolfram kernel is not available in the build environment**, so no *actual* fresh-kernel CG export or B/C coefficient-level ratio is claimed here. The user must run the commands above on the project machine and inspect `comparison.json`. An unsupported CG encoding, source change or failed Matchete stage causes an explicit failure; no synthetic tensors are substituted in B/C reports.

This uses the complex-scalar source Lagrangian directly. It is not yet the independently normalized, fully symmetric real-scalar `lambda_abcd` and Weyl `C_ijab` of Aebischer–Bresciani–Selimović Eq. (4.85). In particular, even if the coefficient ratios are universal, it **does not certify** `GammaCC` normalization, any UV/IR pole assignment, the EFT1–EFT2 counterterm difference, or physical `C5,Y`.

Next step after real data: check the `comparison.json` ratios and implement the independent real-scalar/Weyl conversion with proper operator multiplicities, then compute the renormalization difference rather than assigning a counterterm by negating the raw pole.
