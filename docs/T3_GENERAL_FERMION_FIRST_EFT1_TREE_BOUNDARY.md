# T3 B_-1 / C_-1: fermion-first EFT1 tree boundary

## Scope / verified facts

- Prior **live** user logs: both models passed the native Matchete heavy mass-eigenfield check and exploratory all-heavy one-loop coefficient extraction. Frozen-reference scalar-aligned one-scale masses and SM+Weinberg endpoint running also passed on the user's Windows machine.
- This **new** experiment has **not** yet been run in Matchete. Offline Python tests are not a physical EFT matching certificate.
- Source checked: `Lagrangian/model/T3Fields.wl` supports `Mass -> {Heavy,mass}` / `{Light,mass}`; `Lagrangian/RunMatching.wl` exposes `RunT3Matching[lag,eftOrder,loopOrder]` and implements `Match` plus Greens/EOM/loop/effective-coupling stages.

## Physical threshold intended

The first boundary in the documented fermion-first hierarchy is

    UV(F, X, Y) -- integrate out F at mu_F --> EFT1(X, Y, LLSS)

where `NewFermion` is registered with its original heavy mass `MF`, and both scalar mass eigenstates `NewScalarX` and `NewScalarY` are registered as **Light**, with masses `MSX, MSY`. Their mixed-doublet mass matrix was already diagonalised upstream; the SU(2)-invariant source Yukawa and scalar vertices are reconstructed using native Matchete CG tensors.

The new Wolfram script builds the same 11 ordinary source operators plus 7 additional quartic source operators and verifies `CheckLagrangian`, then runs:

```wolfram
RunT3Matching[lag, 5, 0]
```

The intended tree-level EFT1 Wilson sector has two lepton fields and two **surviving** scalars, with schematic structure `LL Xdag Xdag`, `LL Xdag Ydag`, `LL Ydag Ydag` (and Hermitian conjugates), modulo Matchete's exact gamma/field conventions. No scalar is integrated out at this stage.

## Implemented fail-closed contract

- Depends on the **prior live-certified** native eigenfield result file for the same model, same numerical mass diagonalisation and same benchmark parameter export.
- Requires successful Matchete `CheckLagrangian`, including `DetailedOutput["HeavyMassBasis"]`, with only `F` declared Heavy in the new runner.
- Requires `RunT3Matching[lag,5,0]` to return `Status==Success` and a `MatchedEFT` expression.
- Exports the full matched EFT and selects a diagnostic `GammaCC` + `Field[NewScalarX|NewScalarY]` sector, following the existing `RunMatching.wl` Wilson-seed text convention. Requires at least one term and evidence of both eigenfields, with `NewFermion` field absent.
- Retains the full EFT export even if the diagnostic fails. JSON, progress file and Wolfram logs identify the failure stage; no fabricated Wilson coefficient or matching certificate is written.
- No files from earlier successful pipeline batches are modified or overwritten.

## Files in this package

- `Lagrangian/RunGeneralT3FermionFirstEFT1Tree.wl` (new full Wolfram source)
- `scripts/RunGeneralT3FermionFirstEFT1Tree.py` (new full CLI)
- `tests/python/test_general_t3_fermion_first_eft1_tree.py` (new offline tests)
- `docs/T3_GENERAL_FERMION_FIRST_EFT1_TREE_BOUNDARY.md` (this file)

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_General_Fermion_First_EFT1_Tree_Probe.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_fermion_first_eft1_tree.py
python -m scripts.RunGeneralT3FermionFirstEFT1Tree --model B_-1
python -m scripts.RunGeneralT3FermionFirstEFT1Tree --model C_-1
```

Saved output:

`output/general_t3_mass_eigen_matchete/eft1_fermion_first_tree/`

- `B_m1_final_eft1_fermion_first_tree.json`, `*.progress.json`, `*.log`
- `B_m1_final_eft1_fermion_first_tree_eft1_tree.inputform.txt`
- `B_m1_final_eft1_fermion_first_tree_llss_tree.inputform.txt`
- C analogues; `summary_final.json` after both succeed.

## What remains after this stage

1. Verify first-threshold one-loop hard effects and extract the full EFT1 Wilson set under fixed one-loop power counting; tree `LLSS` alone is **not** final `C5`.
2. Map the native Matchete quartics/Wilson tensor components to the RGBeta/EFT1 conventions with verified normalization; existing source counts do **not** certify physical component maps.
3. Evolve EFT1 from `mu_F` to the scalar thresholds; implement correct mass-eigenstate scalar threshold matching (including distinct mass thresholds if needed).
4. Construct symmetric final `C5`, run SM+Weinberg EFT and compare low-energy masses consistently.
5. Only then extend the matched general-model calculation beyond the two nontrivial mixed-scalar cases.

The previously reported ~0.041 / 0.051 eV low-scale numbers are **one-scale SM-only-background** endpoint results and are not replaced or reclassified by this tree-stage probe.
