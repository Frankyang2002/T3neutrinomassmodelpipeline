# General T3: RGBeta quartic projection fix

## Observed local failure

B_-1 UV: `BetaTerm failed for lambdaT3GeneralQ110101I1` (exit 8).
This coupling shares its field list `{H, Bar[H], Bar[S1], Bar[S2]}`
with the existing `lambdaHHdagS1barS2barCross`, but uses a second candidate
SU(2) contraction.

RGBeta's `Package/Betas.m` implements `CheckPotentialCouplingMixing` and warns
that calling `BetaTerm[oneCoupling, 1]` when two couplings share field content
can return an unpredictable mixture. Its supported alternative is
`BetaTerm[Quartic, 1]`, which uses `ProjectionToUnmixedBetas` to invert the
quartic projection matrix.

## Code change

`RGE/running/rgbeta/RunGeneralT3RGBeta.wl` now calls `BetaTerm[Quartic,1]`
**once** and assigns its *projected* entries to **all** existing quartic
keys as well as new generated couplings. It checks for missing or unevaluated
betas, refuses duplicate names, and records
`GeneralScalarQuarticBetaBasis=RGBetaFullQuarticProjection` in JSON metadata.
Gauge, Yukawa and scalar mass beta entries retain their previous paths.

This avoids a logically inconsistent calculation that would use an unmixed new
quartic alongside mixed old quartic betas. The RGBeta projection can still
fail if the generated contractions are redundant in RGBeta's actual tensor
basis. That must be diagnosed on the local Wolfram installation; the Python
tests alone do NOT validate physical independence, exact normalisations, or
native gauge invariance.

## Test and run

```powershell
python -m pytest -q tests\python\test_general_t3_general_quartic_r_g_beta.py tests\python\test_general_t3_general_quartic_quadratics.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3GeneralQuarticRGBeta.ps1 -ExportBetas -CheckGauge -Models B_-1,C_-1
```

If RGBeta reports `ProjectionToUnmixedBetas::singular`, treat it as a genuine
possible redundant/invalid tensor-basis issue; don't mask it by reverting to
`BetaTerm[symbol,1]`.
