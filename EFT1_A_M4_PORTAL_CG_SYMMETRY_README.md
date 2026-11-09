# A(-4) native portal CG symmetry check

This is a **read-only, additive** diagnostic. It does not alter UV matching, the EFT1 RGE, or the physical Wilson coefficient. It requires the general T3 Lagrangian builder already present in the local project and the local Wolfram/Matchete installation.

## Motivation

Uploaded `RawEFT_lambdaH2Inv1_only.inputform.txt` and `RawEFT_lambdaH2Inv2_only.inputform.txt` each have four equal-weight terms. Relabeling dummy indices in the first expression yields a tensor contribution proportional to `(T_ijab + T_ijba)/2`, with conjugated and unconjugated scalar fields held fixed. If the first CG tensor is antisymmetric in its final two slots, this is exactly zero. That does **not** imply the unsymmetrized Higgs-triplet adjoint portal `H^dagger tau^A H S_2^dagger T^A S_2` is zero.

The prior `GreensSimplify` test demonstrated the isolated `Inv1` contribution simplifies to zero. This new diagnostic tests the native CG tensor symmetry directly and provides independent complex-field contraction probes to distinguish an identically zero CG from a nonzero invariant accidentally annihilated by symmetrization or a phase convention.

## Installation

Extract the ZIP into the root of `C:\thesistemp\T3neutrinomassmodelpipeline`.

## PowerShell commands

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$env:T3_EXPERIMENTAL_GENERAL_FERMIONS = '1'
$env:T3_EXPERIMENTAL_GENERAL_QUARTICS = '1'
$env:T3_EXPERIMENTAL_GENERAL_QUADRATICS = '1'
$env:T3_EXPERIMENTAL_FLAVOR_MASSES = '1'

wolframscript -file Lagrangian\CheckAminus4PortalCGSymmetry.wl `
  output\eft1_general_lagrangian_trace\A_m4_cg_symmetry_01

Get-Content 'output\eft1_general_lagrangian_trace\A_m4_cg_symmetry_01\portal_cg_symmetry.json' -Raw
```

## Expected interpretation

- `inv1.purely_antisymmetric_scalar_slot_pair = true` plus a nonzero `inv1.max_abs_tensor` shows that the 4-term symmetric projection eliminates a **nonzero** native CG invariant. It does **not** prove that Matchete alone introduced this symmetrization; the `PlusHc` construction and original source need checking before deciding on a fix.
- `inv1.max_abs_hermiticity_difference` and `independent_complex_field_probes[*].inv1_plus_hc_for_real_coupling` help diagnose whether the original CG requires an explicit phase for a real hermitian coupling.
- `inv2` is an independent control.

Do **not** add the first portal back as a duplicate term or remove the general-Lagrangian gate before the index orientation and normalization are certified. The current numerical A(-4) EFT1 RGE already evolves 17 RGBeta couplings, including both `lambdaH2` and `lambdaH2Adj`, but this fact does not establish the Matchete/RGBeta basis dictionary.

### Project paths

- `Lagrangian/CheckAminus4PortalCGSymmetry.wl` (new)
- `EFT1_A_M4_PORTAL_CG_SYMMETRY_README.md` (new)

No existing files are replaced.
