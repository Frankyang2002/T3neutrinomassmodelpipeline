# T3 EFT1 full seven-term + one simultaneous two-scalar matching: integrated bridge

## Status, dependencies, and exact project paths

New files only; **no existing scripts, numerical outputs or physics source files replaced**:

- `RGE/running/intermediate/IntegratedEFT1JointScalarBridge.py`
- `RGE/running/intermediate/README_IntegratedEFT1JointScalarBridge.md`
- `scripts/RunIntegratedEFT1JointScalarBridge.py`
- `tests/test_integrated_eft1_joint_scalar_bridge.py`

The code consumes the already installed `SourceBackedEFT1Batch`, `ProvisionalSourceKernel`, `GeneralEFT1C5Mixing`, `FullSU2SourceResponse`, and the **existing** `Numerical/joint_scalar_matching/core.py`, `large_batch.py`, `fermion_boundary_bridge.py`. The joint scalar matcher is **imported**, never copied or independently implemented in the shipped code.

This package was prepared after inspecting read-only GitHub `main` on 2 Nov 2026, SHA `3bb2c659e422cb007f6bf362be3a7dc8e3e5f350`. The repository's source release and Stage 22 extra-quartic geometry were not physically certified at that commit. The new bridge cannot supply missing physical couplings or the Matchete `GammaCC` factor. All numerical outputs are explicitly **per unit unknown signed spinor normalization N_Gamma**.

Fixed matching sequence, with **one** scalar threshold:

`UV -> integrate out ALL fermions at mu_F -> EFT1(H,S1,S2) -> integrate out ALL scalar components together at mu_S -> SMEFT`.

## What is implemented

1. Check all 16 published F-tree LLSS release manifest entries against the exact source file and SHA-256. Preserve every PL source, including numerically cancelling source channels. No PR conjugates are added again.
2. Source-to-master tensor bridge, full complex flavour, with the same source/matching orthogonal real-scalar basis transformation `s_source = U s_matching` used by the existing joint scalar matcher.
3. Strict equality (for this fixed-order diagnostic) of **every complete HHSS component and each named independent quartic channel** in the master and joint scalar inputs, after the explicit basis transformation. Different input scales, dropped portals, duplicate names or mixed-basis errors fail before writing.
4. Compute all **seven** master EFT1 one-loop source terms; demand **full** SU(2) Weinberg-subspace closure of their *sum* and of the quartic part, not only the neutral pivot.
5. Reuse the repository's simultaneous joint bubble with complete real heavy-scalar mass mixing and unequal/degenerate eigenstate masses. Demand SU(2) closure of the finite bubble.
6. Add fixed-order EFT1 source log plus the joint bubble **exactly once** and report the increment in the *master* convention per unknown `N_Gamma`. Include separate matrices for all PL sources, quartic channels, heavy-scalar eigenstate pairs, plus a 21-point source-log trajectory in the upstream seven-term diagnostic (not copied to this compact final output).
7. Independently check `mu_S` cancellation of the scalar-quartic log and the scalar bubble at frozen couplings and masses. **Do not** assert cancellation of the other master terms: their matching-scale partition involves light-field anomalous dimensions and the running of other EFT parameters.
8. Without numerical input, generate 16 source-checked *contracts* (the command that works immediately). No physical or toy quartic is manufactured.

## Factor of two: check rather than assume

The project's master RGE uses

\[
  {\cal L}_5 = \tfrac12 C_5^{pq}(L_p\cdot H)(L_q\cdot H)+\mathrm{h.c.},\qquad
  T_{\mathrm{master},ijAB}=\tfrac14(z_{iA}z_{jB}+z_{iB}z_{jA}),
\]

whereas the current `Numerical/joint_scalar_matching/core.py` uses a Weinberg template

\[
 T_{\mathrm{joint},ijAB}=\tfrac12(z_{iA}z_{jB}+z_{iB}z_{jA})=2T_{\mathrm{master},ijAB}.
\]

**Therefore `C5_master = 2*C5_joint` for the same external LLHH tensor.** The bridge projects each scalar bubble directly with the master template and independently verifies equality to twice the existing joint-matcher coefficient. No duplicate scalar-bubble amplitude is added.

The fixed-order source-induced increment, not including independent direct hard-F matching, is

\[
  \left.\Delta C_5^{\mathrm{SMEFT}}(\mu_S)\right|_{\mathrm{LLSS}}
 = N_\Gamma\left[\frac{\ln(\mu_S/\mu_F)}{16\pi^2}
  \beta_{C_5}^{(1),7\mathrm{-term},\mathrm{formal}}
 + \Delta C_{5,\mathrm{joint\ bubble}}^{\mathrm{formal,master}}(\mu_S)\right].
\]

The unknown `N_Gamma` remains JSON `null`; it is **not silently set to one physically**. The result is per unit unknown normalization, not a certified mass matrix. `C5(mu_F)` is **not supplied or assumed zero**. Complete hard heavy-fermion one-loop matching, possible other dimension-five sources, independent full one-loop scalar threshold counterterms, and SM running beyond the threshold are excluded.

The bubble uses O(hbar^0) LLSS inserted into one O(hbar) loop; LLSS self-running inside that loop would be O(hbar^2) and is not double counted.

## PowerShell: ready now, all 16 input contracts

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive -LiteralPath "$HOME\Downloads\T3_EFT1_SevenTerm_JointScalar_Batch.zip" -DestinationPath .
python -m unittest discover -s tests -p 'test_integrated_eft1_joint_scalar_bridge.py' -v

python scripts\RunIntegratedEFT1JointScalarBridge.py `
  --release-root 'output\general_observable\lo_fermion_boundary\releases\source_021dcc1e8fe8ba0f' `
  --output-root 'output\general_observable\lo_eft1_mixing\integrated_eft1_joint_v1'

Get-Content 'output\general_observable\lo_eft1_mixing\integrated_eft1_joint_v1\summary.json' -Raw
```

Expected: 16 `integrated_input_contract.json` files under model subdirectories, and `summary.json` with `model_count=16`, `models_evaluated_conditionally=0`, `physical_C5_computed=false`.

Run with `--check-only` first if you want to validate the inputs without any writes. **Never overwrite existing output directories**: use a distinct `--output-root` if a run already exists.

## PowerShell: conditional numerical integration after external inputs exist

Supply, **read-only**, the same model-matched inputs required by the existing two workstreams:

- `<coupling-root>\A_m4\couplings.json` etc: `t3-eft1-explicit-master-couplings-diagnostic-v1`, using complete `Y`, all independent `V4` quartic channels, and SU(2)/U(1) gauge tensors (or explicitly toy-zero diagnostics). The source SHA-256 and renormalization scales must match. Set the explicit `coverage_assertion.fixed_order_quartic_identical_muF_muS=true` **only if that numerical assumption is justified** for this frozen fixed-order run.
- `<background-root>\A_m4.json` alongside `<background-root>\A_m4.npz` etc: `t3-joint-scalar-conditional-background-v1`, with real orthogonal `U`, a positive **full** real mass-squared matrix of both scalar multiplets at one `mu_S`, and the full `lambda_HHSS_muS` + independently labelled quartic channels. `source_scalar_real_labels` must match the actual F release, including X/Y and split-real cases.

Every named quartic channel must be present in **both** inputs; a representative/candidate from Stage 22 is **not** a physical coefficient. The bridge forbids EWSB-resummed masses at this dimension-five Weinberg-check step. If full source-specific numerical couplings are missing, this numeric mode blocks rather than inserting zeros.

```powershell
python scripts\RunIntegratedEFT1JointScalarBridge.py `
  --release-root 'output\general_observable\lo_fermion_boundary\releases\source_021dcc1e8fe8ba0f' `
  --coupling-root 'output\general_observable\lo_eft1_mixing\external_diagnostic_couplings' `
  --background-root 'output\general_observable\lo_joint_scalar_matching\scalar_backgrounds' `
  --output-root 'output\general_observable\lo_eft1_mixing\integrated_eft1_joint_with_backgrounds_v1'
```

For a single model add `--model B_-1`. All requested models are validated **before** any results are written; failing physical input means no partial published output. Every conditional result is `output/general_observable/lo_eft1_mixing/<new-run-name>/<tag>/conditional_eft1_joint_scalar.json`. `physical_C5_computed` and `physical_complete_one_loop_certified` remain false in all such results.

## Test scope

The new tests cover SU(2) projection against a deliberately constructed **synthetic** complex Weinberg tensor, real basis rotations, independent heavy-scalar channels, exact Weinberg-operator factor-of-two, nonzero one-loop log and scalar bubble, matching-scale-log cancellation, 16-model safe preflight and refusal to overwrite, missing/faulty quartet background, missing quartic channel and missing frozen-coupling declaration. In environments where the joint scalar workstream is not mounted, tests use an explicit local algebraic stand-in solely to test the bridge; they do **not** present that stand-in as verification of the repository implementation. The installed project's tests use its real joint-scalar code.

**Next genuine physics blocker:** finish the independently sourced all-16 EFT1 coupling dictionary (rather than inserting Stage 22 candidate symbolic geometry), and separately certify the Matchete signed spinor amplitude. Until then the immediately runnable mode is an input audit, not a neutrino-mass calculation.
