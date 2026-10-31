# Workstream C — original fermion-hard LLHH boundary (v13.1: import correction)

## Why this matters

The v12 interface provides a conditional, dimension-five **Y scalar-bubble** contribution but does not include the hard one-loop contribution generated when integrating out the fermion F. The original repository already contains those F-first native expressions:

- `output/general_t3_mass_eigen_matchete/eft1_fermion_first_oneloop/B_m1_final_eft1_fermion_first_oneloop_raw_hard_llhh_sector.inputform.txt`
- `output/general_t3_mass_eigen_matchete/eft1_fermion_first_oneloop/C_m1_final_eft1_fermion_first_oneloop_raw_hard_llhh_sector.inputform.txt`

The original saved GitHub blobs were read-only inspected before implementation. **Both** contain 204 additive terms: 68 pole, 68 `Log[μbar2/MF^2]`, and 68 scale-independent terms. Each sector has 34 PR and 34 PL terms. This is the first concrete missing input to the sequential C5 handoff. The user should use the **current local** files; the audit records their hashes and refuses changed inventories rather than silently replacing them with GitHub copies.

## What the script computes

The script reuses the current strict native `Numerical.experimental.T3BCNativeYLLHHFlavorKernel._record_from_term` parser to extract exact neutral SU(2) scalar/flavour/spurion coefficients of all three **F-hard source sectors**. It reuses v12's ordered flavour/Warsaw normalization, gated on v11.2's installed Matchete `ExactPlusEquality`. It checks:

1. Exact pole = `Log[μbar2/MF^2]` source coefficient, first as opaque InputForm monomials and then as projected neutral components.
2. Native source chirality, Higgs and epsilon orientation, MF and hbar powers, Yukawa slots, and permitted quartics.
3. Exact PR/PL Hermitian-conjugate coefficients.
4. Correct `A_PL -> C5 = A_PL + A_PL.T` flavour factor without counting PR twice.
5. Whether the *observed* F constant polynomial equals the F pole polynomial; this is a diagnostic only, not a counterterm assertion.

The output preserves **three separate F-boundary sectors**. It references v12 as a different threshold, without adding these expressions or treating them as physical Wilson coefficients.

## Physics interpretation and limits

At the F boundary, raw one-loop `LLHH` terms contain coefficients associated with a heavy-F expansion. A formal finite term at `μF=MF` is **not** a renormalized matching coefficient until the UV/EFT1 amplitude subtraction, UV/IR separation, wavefunction/operator counterterms and scheme are fixed. The same one-loop UV diagram may have distinct hard F, Y and X momentum regions; the source terms must **not** simply be numerically added at an arbitrary common scale.

The output is separate from Y, from inherited `LLXX`/`LLHX`, from EFT1 running, and from X matching. No `physical_C5`, neutrino mass, or independent F counterterm is certified.

The v12 scalar-bubble script and original matching/Lagrangian/RGBeta code are not modified.

## Project installation

The ZIP contains the complete script at:

- `scripts/BuildT3BCFHardBoundaryInterface.experimental.py`

and this README at:

- `docs/WorkstreamC_FHardBoundary_Interface_README.md`

From the repository root:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WorkstreamC_FHardBoundary_v13.zip" -DestinationPath . -Force

python scripts\BuildT3BCFHardBoundaryInterface.experimental.py --self-test
python scripts\BuildT3BCFHardBoundaryInterface.experimental.py `
  --output output\audits\t3_bc_f_hard_boundary_v1\summary.json
```

Required earlier local evidence:

- `output/audits/t3_bc_y_conditional_finite_interface_v1/summary.json` (v12)
- `output/audits/t3_bc_y_warsaw_direct_identity_v11/summary_live.json` (v11.2)

The script also reads original F-first B/C native `raw_hard_llhh_sector.inputform.txt` files at the paths above, plus the existing v12 script and native kernel parser.

## Next *physical* operation

Compute an independently renormalized F-to-EFT1 matching difference, including the EFT1 soft subtraction and inherited coefficient's scale evolution. Then combine with the Y matching in a single consistently normalized and scheme-correct EFT sequence. The present output identifies and exports the missing F-hard boundary data but **does not solve this subtraction**.

## v13.1 correction (October 2026)

v13's local synthetic self-test passed but the actual run stopped with
`No module named 'Numerical'`. Executing a file in `scripts/` does not
necessarily add the repository root (containing `Numerical/`) to
`sys.path`.

This full replacement calls `check_project_imports(root)` in **both**
`--self-test` and the actual run. It checks that the two required native
modules exist in `Numerical/experimental/`, prioritizes the explicitly
selected local root, and verifies that the modules are loaded from those
exact files rather than unrelated site-packages.

No coefficient extraction, loop subtraction, source-file selection or
physical-certification logic has been changed. The self-test now includes
the import path used by production.

Install `T3_WorkstreamC_FHardBoundary_v13_1.zip` over the repo, then run:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WorkstreamC_FHardBoundary_v13_1.zip" -DestinationPath . -Force
python scripts\BuildT3BCFHardBoundaryInterface.experimental.py --self-test
python scripts\BuildT3BCFHardBoundaryInterface.experimental.py `
  --output output\audits\t3_bc_f_hard_boundary_v1\summary.json
```

The extra `PYTHONPATH` workaround is no longer required.
