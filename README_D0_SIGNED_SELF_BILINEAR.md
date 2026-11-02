# D(+0) shared-Yukawa sensitivity: signed-real self-bilinear mass fix

## Scope

One full replacement: `Numerical/matching/GeneralT3FlavourC5.py`.
One new offline test: `tests/python/test_thesis_signed_self_bilinear.py`.
No edits to the Wolfram/Matchete source, general RGEs, matching coefficients, or shared A(-2) report generator.

The log reports that D(+0)'s self-bilinear rejects the evaluation at gamma >= 1,
while some points below gamma=1 evaluate. It does not print the value of the
running `muT3GeneralS1Sq`; a sign reversal is a plausible, not proven, cause.
The updated guard reports the actual complex evolved value if it rejects it.

### Fixed Matchete real-field convention

The installed `Lagrangian/RunT3GeneralSelfMassBasis.wl` uses
`S = (S_minus + i S_plus)/sqrt(2)` and matches independent real eigenfields.
With the Lagrangian's holomorphic mass convention,

- `m_minus^2 = m_S^2 - 2 Re(mu^2)`
- `m_plus^2  = m_S^2 + 2 Re(mu^2)`

The existing helper `GeneralT3MultiModelNeutrinos.mass_params` instead sorts
these masses using `abs(mu^2)`. This agrees only when `mu^2 > 0`.

For *negative real* `mu^2`, this replacement swaps the numeric values
assigned to the already-defined `PlusMass` and `MinusMass` Matchete symbols.
It **does not** swap field variables, change vertices, discard a bilinear or
silently treat an arbitrary complex phase as real.

The accepted imaginary part is roundoff-level only, preserving the previous
`1e-10 * max(1,abs(mu^2))` imaginary tolerance. A genuinely complex bilinear
still fails with a diagnostic that includes the actual value.

This is still a provisional one-scale general Matchete C5 matching model,
not certified sequential matching or a physical Matchete/RGBeta normalization.

## Windows PowerShell

From project root:

```powershell
Expand-Archive "$HOME\Downloads\T3_D0_Signed_Real_Bilinear_Fix.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_thesis_signed_self_bilinear.py
python -X utf8 -m scripts.RunT3SharedA2General16Reports `
    --donor-run-id thesis16_general_c5_full_v1 `
    --run-id shared_A2_general16_ovals_v1 `
    --gamma-min 0.7 --gamma-max 1.4 --gamma-points 7 --resume
```

`--resume` reuses successful scans for 15 models and recomputes D(+0).
On success, the existing reporter re-creates `sensitivity_oval.html` with 16
models. If D(+0) fails for a genuine complex phase, read its per-point error:

```powershell
$dir = 'output\studies\general_scalar_extended\shared_a2_general16_reports\shared_A2_general16_ovals_v1\D_p0'
$scan = Get-Content "$dir\sensitivity_diagnostics.json" -Raw | ConvertFrom-Json
$scan.directions.T1.points | Select-Object gamma,status,message | Format-Table -Wrap
```

Do not replace negative `mu^2` by its absolute value without the correct
field-label mapping, and do not discard imaginary parts exceeding the roundoff
tolerance. Full complex-phase matching requires a separate derived rotation.
