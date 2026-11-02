# Conditional RGBeta/Matchete portal-sign sensitivity

The all-32 native tensor identity audit found that the RGBeta portal tensor
`delta_H*delta_S` is in the real linear span of the Matchete portal CGs.
This constrains tensor geometry, but not the physical cross-package coupling map.

The source Matchete portal Lagrangian contains `+1/2 PlusHc[lambdaHInvk*T_k]`.
The RGBeta model-builder convention defines the quartic part of the
Lagrangian with a minus sign (Thomsen, Eur. Phys. J. C 81, 408 (2021),
Secs. 2, 4; https://doi.org/10.1140/epjc/s10052-021-09142-4).
Consequently, for a common real field and CG basis, a **candidate** physical
identification is `lambdaHInvk_Matchete = -c_k*lambdaH_RGBeta`, where `c_k`
are the measured all-four-index identity-span coefficients.

The existing **uncertified** direct UV adapter uses
`A_ii=mSiSq-lambdaHi*v^2/2` (Matchete-positive-sign interpretation).
This script independently evaluates an RGBeta-sign candidate
`A_ii=mSiSq+lambdaHi*v^2/2`, changing **only** that portal diagonal.
It does not resolve the sign of the lambdaT3 scalar mixing, Yukawa phases,
fermion mass normalization, or symmetry-allowed extra interactions.
Both outputs are provisional. It does not edit prior files.

Install paths:

- `scripts/RunT3DirectOneLoopPortalSignSensitivity.py`
- `tests/python/test_t3_direct_one_loop_portal_sign_sensitivity.py`
- `Numerical/direct_one_loop/PORTAL_SIGN_SENSITIVITY.md`

Run from repository root (PowerShell):

```powershell
$base = 'output/general_observable/lo_direct_one_loop'
python -m unittest discover -s tests/python -p 'test_t3_direct_one_loop*.py' -v
python scripts/RunT3DirectOneLoopPortalSignSensitivity.py `
  --project-root . `
  --input-dir "$base/provisional_inputs" `
  --baseline-dir "$base/provisional" `
  --identity-audit "$base/source_audits/native_portal_identity_decomposition.json" `
  --output-dir "$base/portal_sign_sensitivity"
Import-Csv "$base/portal_sign_sensitivity/all_16_portal_sign_sensitivity.csv" `
  | Select-Object model, baseline_m3_ev, variant_m3_ev, relative_m3_change `
  | Format-Table -AutoSize
```

Output is restricted to
`output/general_observable/lo_direct_one_loop/portal_sign_sensitivity/`.
No EFT, matching, RGE, threshold, frozen benchmark, earlier CG audit or
previous one-loop input/output is modified. Numerical change is a sensitivity
probe, **not** a certified correction.
