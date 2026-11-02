# Batch 1: native split-triplet Ward repair and verification

**Status:** all four original v6 Weinberg projection failures resolved in the
existing corrected all-16 GitHub numerical summary. This bundle **repackages,
without modification**, the two full replacement files already present in the
latest read-only GitHub repository. No GitHub write was made.

## Evidence

- Matchete source: `Lagrangian/RunT3GeneralSelfMassBasis.wl` defines
  `NewScalar[k] = (T3ScalarMinus[k] + i T3ScalarPlus[k])/sqrt(2)` with unchanged
  native SU(2) index `k` (no spherical intertwiner).
- Replacement source:
  `Numerical/joint_scalar_matching/special_mass_v6.py`, Git blob SHA-1
  `2fb97f9eda33b1780336124f5380d085b51ab240`.
- Replacement regression tests:
  `tests/test_t3_joint_scalar_special_mass_v6.py`, Git blob SHA-1
  `2aae0eadc9cfe4a0aca90b4680f13497b9900274`.
- Existing corrected result:
  `output/general_observable/lo_joint_scalar_matching/exploratory_all16_v6_native_fix/summary.json`,
  Git blob SHA-1 `53817161e7fee40562a8834b437d8217407d92cd`.
- Current corrected summary: `models_checked=16`,
  `models_with_conditional_C5=16`, `errors=0`.
- Original source bug: spherical prototype `W` was applied to native Matchete
  Cartesian real components of the **triplet**. The same error did not affect
  self-split singlet models. The repair uses an index-preserving permutation
  for real/imaginary scalar components, and builds the full fixed-basis Hessian
  for complex `mu` rather than rerotating the already-matched LLSS coefficient.

## Formerly failing models: corrected numerical evidence

| model | full Weinberg projection residual | C5-related neutrino-mass singular values, eV per abs(N_GammaCC) |
|---|---:|---|
| A_-2 | 3.119401668e-16 | 0.02406817066, 0.02005680867, 0.01805107442 |
| D_+0 | 1.543185790e-16 | 0.02442053435, 0.02035044508, 0.01831534641 |
| E_-2 | 3.099511452e-16 | 0.03960790676, 0.03300658861, 0.02970584192 |
| E_+0 | 1.975106730e-16 | 0.03952674136, 0.03293895078, 0.02964496805 |

These values are taken from the corrected GitHub summary, **not** from an
independent execution of the full project in this packaging container.
The local files in this bundle were checked for Python syntax and exact
Git blob equality. The GitHub documentation reports 58/58 local tests; the
independent full run remains to be verified on the user's Windows checkout.

## Install and run (PowerShell)

From `C:\thesistemp\T3neutrinomassmodelpipeline`, with the downloaded
`T3_Batch1_NativeWard_FullReplacement.zip` placed in Downloads:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Join-Path $HOME 'Downloads\T3_Batch1_NativeWard_FullReplacement.zip'
if (-not (Test-Path -LiteralPath $zip)) { throw "Missing $zip" }
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
python -m unittest discover -s tests -p 'test_t3_joint_scalar_*py' -v
python -m scripts.RunT3JointScalarAll16V6 `
  --allow-exploratory-map `
  --allow-exploratory-special-bases `
  --output-dir 'output\general_observable\lo_joint_scalar_matching\exploratory_all16_v6_native_batch1_local'
```

The output directory is intentionally new; the runner refuses to overwrite
previously validated outputs. Check `summary.json` for 16 conditional C5
results and zero errors. If the local source inputs are different, retain
those diagnostics and do not suppress Ward failures.

## EFT and physical limitations

`UV -> integrate ALL heavy BSM fermions at muF=100000 GeV -> EFT1 with BOTH
BSM scalars -> integrate BOTH scalar multiplets at muS=1000 GeV -> SMEFT`.
SM fermions remain active throughout both BSM thresholds.

**EXPLORATORY / CONDITIONAL — NOT A COMPLETE PHYSICAL ONE-LOOP NEUTRINO-MASS PREDICTION**

The overall signed `N_GammaCC` normalization is unresolved; native
Matchete/RGBeta quartic aliases, the muF/muS scalar-basis alignment, and
hard-fermion matching/full EFT1 anomalous-dimension contributions are not
certified. Scalar masses and couplings are illustrative EFT1 benchmarks.
Neutrino-mass convention: `mnu=(246.22 GeV)^2 C5 / 2`, reported per
`abs(N_GammaCC)`; no PMNS fit is claimed.
