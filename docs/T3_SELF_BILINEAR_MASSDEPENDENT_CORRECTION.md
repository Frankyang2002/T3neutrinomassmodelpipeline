# Corrected six-model A/D/E self-bilinear one-scale matching

## Why the earlier `self_anomalous_v1` contribution was identically zero

Its `unit_topology_kernels` used `scripts.RunT3GeneralScalar48Campaign.matching`, which **prefers** `output/full/hypercharge/T3_<class>_alpha_<tag>/data/final_weinberg_coefficient.json` (the fixed-order hierarchical matching artifact) over the one-scale legacy coefficient.

For A(+0), A(-2), D(+0), D(-2), E(+0), E(-2), all six saved hierarchical expressions have **no MS1 or MS2 dependence**. Substituting `m_plus` or `m_minus` in that kernel returns the same result exactly. Thus `C5_anomalous = coupling * (U_plus - U_minus)/2 = 0` is a **kernel limitation**, not a physical zero. The six Stage29 initial HH extra quartics were nonzero, and the user initialized the quadratic bilinear to 1% of the diagonal mass squared.

## Correction

This new runner uses each model's saved **mass-dependent** original Matchete one-scale `c5_coefficient.txt` through the existing `LegacyKernel` class. It requires explicit `MS1` and `MS2` source symbols and a **nonzero mass response** for nonzero mass splitting. It evaluates the source topology with `lambdaT3=1` at `m_plus` and `m_minus`, then constructs:

- `C5_normal = lambdaT3 * (U_plus + U_minus)/2`;
- `C5_anomalous = conj(lambda_extra_M) * P_alpha(mu2) * kappa * (U_plus - U_minus)/2`;
- `C5_total = C5_normal + C5_anomalous`;

where `lambda_extra_M = -lambda_extra_RGBeta` is **only an assumed sign map**, `P_alpha` is the phase from the bilinear, and `kappa` is checked from saved Matchete native CG tensors. The existing UV-evolved heavy Yukawas/MF, EFT1 endpoint and 1 TeV -> 100 GeV SM+Weinberg running are reused.

**Important:** this changes the matching approximation compared with the old hierarchical calculation. The one-scale source includes scalar masses, but is **not** the full sequential Wilson matching. Its numerical neutrino masses are conditional. The new mixed 16-row table cannot be treated as predictions in a common certified renormalisation scheme. No physical Matchete–RGBeta overall coupling map is claimed.

## Files and installation

Extract into `C:\thesistemp\T3neutrinomassmodelpipeline` (preserving directory structure):

- `Numerical/experimental/T3SelfBilinearMassDependentKernel.py` — NEW
- `scripts/RunT3SelfBilinearMassDependentC5.py` — NEW
- `tests/test_t3_self_bilinear_massdependent_c5.py` — NEW
- `docs/T3_SELF_BILINEAR_MASSDEPENDENT_CORRECTION.md` — NEW

None of the previous scripts/reports or root README are overwritten.

## PowerShell

```powershell
python -X utf8 -m pytest -q tests/test_t3_self_bilinear_massdependent_c5.py

python -m scripts.RunT3SelfBilinearMassDependentC5 `
    --all `
    --run-id continuous_v1
```

To resume **without changing the implementation/input hashes**:

```powershell
python -m scripts.RunT3SelfBilinearMassDependentC5 `
    --all `
    --run-id continuous_v1 `
    --resume
```

If all six outputs are present and only a combined table is needed:

```powershell
python -m scripts.RunT3SelfBilinearMassDependentC5 `
    --all `
    --run-id continuous_v1 `
    --combine-only
```

Output location:

`output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/self_anomalous_massdep_v1/`

The new individual files are `A_p0_massdependent_c5.json`, etc., and the mixed-source table is `combined_16_summary.csv` with explicit matching-category labels. Do **not** reuse the earlier six zero-anomalous values as physical constraints.

## Validation scope

The tests cover mass-source selection, zero-coupling and complex-phase limits, actual scalar-mass sensitivity, and rejection of mass-insensitive or nonlinear kernels. Local isolated tests cannot evaluate the full project Matchete source through the complete user's Windows installation; run the commands above and inspect the printed six ratios. The result is still conditional on a cross-package sign/phase convention and the one-scale approximation.
