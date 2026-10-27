# T3: Matchete scalar operator audit with installed-version baseline calibration

## Source of the latest failure

The 43/43 Python tests from the previous ZIP passed, but the user's **live**
Matchete `CheckLagrangian[..., DetailedOutput -> True]` returned:

- ordinary baseline `CheckLagrangian[base]` **True**, established earlier by
  `BuildT3Lagrangian`;
- explicit `Hermiticity`, `GaugeInvariant`, `CanonicallyNormalized`,
  `HeavyMassBasis`, and `AllObjectsDefined` checks **True**;
- `UncontractedIndices`, `HeavyTadpoles` and spin-chain failures **empty**;
- **`Charges`** rather than the `ChargeNeutral` key assumed by the previous parser;
- `GaugeAnomalies` and `ZeroModeMixing` carrying native values other than
  literal `True` or `{}`.

The previous parser mislabelled these last three as `FailedOrUnsupported`,
causing `BaselineDetailedCheckFailed` even in A_+0. This was a parser problem;
no new physics discrepancy was established.

## Corrected acceptance criteria

The replacement `Lagrangian/RunGeneralQuadraticMatcheteAudit.wl` performs two
actual Matchete calls for each candidate Lagrangian:

1. `CheckLagrangian[lag]` (the ordinary authoritative full check), and
2. `CheckLagrangian[lag, DetailedOutput -> True]` (per-check diagnostic data).

For the original **unextended baseline** it requires the ordinary check to be
`True` and explicit structural diagnostics to pass. It saves exact `InputForm`
values of native diagnostic entries whose pass/fail semantics have not been
established from the installed Matchete version.

For each **individual additional quartic**, **individual quadratic** and the
**combined extra scalar Lagrangian** it demands:

- expected ordinary check: `True` unless the extra B/C(-1) scalar mixing is
  present, in which case the ordinary check must be `False`;
- explicit passing values for known structural tests: Hermiticity, no open
  spin chains or uncontracted indices, canonical normalization, no heavy
  tadpoles, all objects defined, gauge invariance;
- `HeavyMassBasis == True` for ordinary cases and `False` for the mixed
  B/C(-1) quadratic cases;
- the **exact same** native opaque `InputForm` diagnostic values as the
  independently passing baseline, including `Charges`, `GaugeAnomalies`
  and `ZeroModeMixing` when exposed by that installation.

A new/unrecognised native diagnostic key, or a different opaque value,
causes a **failure**, not an automatic pass. Individual and combined reports
include both the normal Matchete result and the full recorded native values.

This is deliberately a **baseline-anchored consistency check**, not a proof
that the semantic meaning of every version-specific Matchete return object
has been decoded. In particular the B/C(-1) mixed operator still cannot be
matched in the original non-diagonal heavy mass basis.

## Files in the ZIP

- `Lagrangian/RunGeneralQuadraticMatcheteAudit.wl` — full replacement
- `scripts/RunGeneralT3QuadraticMatcheteAudit.py` — full replacement
- `tests/python/test_general_t3_matchete_mass_basis.py` — full replacement
- `docs/T3_MATCHETE_DETAILED_MASS_BASIS_AUDIT.md` — this full replacement

No production physics files or GitHub code are modified.

## Commands (PowerShell)

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Expand-Archive -LiteralPath "C:\path\to\T3_Matchete_Baseline_Anchored_Audit_Fix.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_matchete_mass_basis.py
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --model A_+0
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --model B_-1
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --model C_-1
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --all
```

Generated diagnostics are written to
`output/general_t3_quadratic_eft1/matchete/`.

## Certification boundary

A successful experimental audit should report 14 ordinary heavy-basis checks
and two (`B_-1`, `C_-1`) baseline-anchored operator checks with
`HeavyMassBasis` specifically blocked. No Matchete-RGBeta physical coefficient
normalization, mixed scalar mass eigenbasis, or threshold matching is
certified. Existing numerical EFT1 evolutions are not modified.

If the native `Charges` / `GaugeAnomalies` / `ZeroModeMixing` entries legitimately
vary between a baseline and an extended Lagrangian, the audit rejects the
comparison and prints both versions as `InputForm`. That difference must be
analysed rather than forcibly whitelisted.
