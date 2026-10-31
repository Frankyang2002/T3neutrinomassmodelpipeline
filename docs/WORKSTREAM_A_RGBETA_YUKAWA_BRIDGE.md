# Workstream A — six generalized RGBeta Yukawas: registration versus Ward certification

This replaces the previous (native-Ward-required) experiment after the observed
RGBeta v1.2.0 diagnostic on `A_m2`.

## Observed installed-package result (31 October 2026)

`AddYukawa[..., CheckInvariance -> True]` emitted `AddYukawa::gaugeInv`
for U1Y, SU2L **and SU3c** on `A_m2`, with unreduced `TfLeft`, `TStructure`
and `sDel` terms. All interacting fields are SU3c singlets; their U1 charges
also sum to zero. Consequently this is **not** a verified physical gauge
violation. Nor is it a verified Ward identity: the native checker failed to
reduce the symbolic expression to zero. Do not patch internal RGBeta
`OptionTest` or silence `CheckInvariance` errors and call it a pass.

## Isolated registration check

The Wolfram companion builds the ordinary five-Yukawa RGBeta model, then adds
exactly one extra coupling with the default `CheckInvariance` option **omitted**:

| Models | Extra coupling | Field list in RGBeta | Neutral scalar dimension |
|---|---|---|---|
| A_m2, D_m2, E_m2 | `y2Conj` | `{Bar[S2], l, Fc}` | 3, 1, 3 |
| A_p0, D_p0, E_p0 | `y1Conj` | `{S1, l, F}` | 1, 3, 3 |

For both cases, the U1 charge sum is exactly `-1/2 + 0 + 1/2 = 0`.
The SU2 singlet multiplicity in `2 x 2 x (1 or 3)` is one. This is
representation-content eligibility, **not** a verification of RGBeta's
native triplet CG component covariance. The old Matchete report independently
accepted the generalized Yukawa in the full model.

The Wolfram runner verifies the RGBeta `$yukawas` association changes from 5
to 6 entries and independently reads the registered `Fields`, `Invariant`,
`Chirality` and flavour `Indices`. The SU2 invariant function is the same as
the baseline T3 Yukawa at a neutral scalar; this equality does **not**
certify a cross-package triplet component map.

No beta functions are evaluated, and no RGE, matching or numerical source is
modified. The tool reports both `native_ward_check_invoked=false` and
`native_ward_check_passed=false` intentionally.

## Execution

Run at repository root, after backing up any existing output:

```powershell
python Lagrangian\RunT3RGBetaExtraYukawaBridge.py --self-test
python Lagrangian\RunT3RGBetaExtraYukawaBridge.py --all --overwrite
```

Output directory:

`Reports/output/workstream_a/rgbeta_extra_yukawa_bridge/`

- `native_rgbeta/A_m2.json` etc.: actual native registration records
- `native_rgbeta/A_m2.log` etc.: raw Wolfram logs
- `rgbeta_extra_yukawa_bridge.json` and `.md`: aggregate conditional result

Target status if six native registrations succeed:

`SixExtraRGBetaYukawasRegisteredIndependentGaugeCensusPassedWardUncertified`

The result **does not** mean gauge Ward identities, physical Yukawa coupling
normalization, RGEs, or the final UV cross-package dictionary are certified.
The next physical bridge would need a direct component-level test of RGBeta's
scalar/fermion representations and Matchete's registered CG tensors.
